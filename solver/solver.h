
/* Solves for approximation of given target kernel. */

#ifndef _solver_h
#define _solver_h

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <set>
#include <map>
#include "array.h"
#include "timer.h"
#include "pareto.h"
#include "filter.h"
#include "params.h"
#include <lbfgs.h>
#include <new>

extern adept::Stack stack;
using std::pair;
using std::set;
using std::map;
using std::dynamic_pointer_cast;
//using namespace std;

#define SPARSE_NULL shared_ptr<FilterSparse<real> >(NULL)

#define ATTEMPT_2MOVES         0
#define CHECK_NAN              0
#define RECOVER_NAN            0
#define PENALIZE_LOW_FSCALE    0
#define LOW_FSCALE_PARAM       1e-4
#define FAST_CONTRACT_RANDOM   1
#define DEBUG_TRANSLATE_BORDER 0
#define DEBUG_MULTIBAND        0

#define OPTIMIZE_SKIP          1e150          /* Error value indicating optimize_and_add() already has seen mask_samples of sparsity mask */

#define RECLAIM_ADEPT_MEMORY() stack.new_recording()

extern vector<vector<int> > antialias_subsample_x;      /* Antialiasing (dx, dy) patterns */
extern vector<vector<int> > antialias_subsample_y;

/* ------------------------------------------------------------------------
   Pareto frontier type for sparse filters: use PARETO_TYPE(real).
   ------------------------------------------------------------------------ */

class Profiling { public:
    double T_fit;
    double T_deriv;
    double T_bfgs;
    int total_of_calls;         /* Total calls since solver start */
    int current_iter;
    int of_calls;               /* Calls within current contract, expand, translate step */
    int bfgs_calls;
    int bfgs_iters;
    int num_solve_initial;
    int num_solve_add;
    int add_count;
    int iter_changed;
    int max_iterations;
    int skip_count;
    map<void *, int> adaptive_iter;   /* Last used on given iteration */
    map<string, int> mask_samples;
    void init();
    void reset();
};

extern Profiling profiling;

/* ------------------------------------------------------------------------
   Fitting functions
   ------------------------------------------------------------------------ */

template<class real>
real kernel_sum(const Array<real> &a) {
    real ans = 0;
    for (int i = 0; i < a.nelems; i++) {
        ans += a.data[i];
    }
    return ans;
}

#define GET_OBJECTIVE_BBOX() \
    int xmin = 0, ymin = 0, xmax = out.width(), ymax = out.height(); \
    if (params.ignore_boundary > 0) { \
        xmin = ymin = params.ignore_boundary; \
        xmax = out.width()  - params.ignore_boundary; \
        ymax = out.height() - params.ignore_boundary; \
    }

template<class real>
real filter_scale(const Array<real> &out, const Array<double> &target) {
    if (params.scale_correct) {
        GET_OBJECTIVE_BBOX();
        real num = 0, denom = 0;
        if (!params.weights.size()) {
            for (int y = ymin; y < ymax; y++) {
                for (int x = xmin; x < xmax; x++) {
                    real out_v = out(y, x);
                    num += target(y, x) * out_v;
                    denom += out_v * out_v;
                }
            }
        } else {
            Array<double> *W = params.weights_array();
            for (int y = ymin; y < ymax; y++) {
                for (int x = xmin; x < xmax; x++) {
                    double w_v = (*W)(y, x);
                    real out_v = w_v * out(y, x);
                    num += w_v * target(y, x) * out_v;
                    denom += out_v * out_v;
                }
            }
        }
        return num / MAX(denom, real(1e-12));
    } else {
        return 1.0/kernel_sum(out);
    }
}

template<class real>
real norm_target(const Array<real> &a) {
    real sum2 = 0;
    for (int i = 0; i < a.nelems; i++) {
        sum2 += a.data[i]*a.data[i];
    }
    return sqrt(sum2);
}

template<class real>
Array<real> normalize_target(const Array<real> &a) {
    if (params.weights.size() == 0) {
        real fscale = (real(1.0) / norm_target(a));
        return a * fscale;
    } else {
        Array<double> *W = params.weights_array();
        ASSERT(W, "expected weights array non-NULL");
        ASSERT(W->nelems == a.nelems, "expected weights and target of same size");
        real sum2 = 0;
        for (int i = 0; i < a.nelems; i++) {
            real current = a.data[i] * W->data[i];
            sum2 += current * current;
        }
        real fscale = (real(1.0) / sqrt(sum2));
        return a * fscale;
    }
}

/* Quantized sqrt() */
double qsqrt(double E);

//void get_spatial_weights(const Array<double> &target, Array<double> *&w_ptr);

bool is_symmetry();

template<class real>
real error_given_fscale(const Array<real> &out, const Array<double> &target, real fscale) {
    real ans = 0;
    GET_OBJECTIVE_BBOX();
    if (!params.weights.size()) {
        for (int y = ymin; y < ymax; y++) {
            for (int x = xmin; x < xmax; x++) {
                real delta = target(y, x) - out(y, x)*fscale;
                ans += delta * delta;
            }
        }
    } else {
        Array<double> *W = params.weights_array();
        ASSERT(W, "expected W to not be NULL");
        ASSERT(W->width() == target.width(), "expected weights width == target width");
        ASSERT(W->height() == target.height(), "expected weights height == target height");
        for (int y = ymin; y < ymax; y++) {
            for (int x = xmin; x < xmax; x++) {
                real delta = (target(y, x) - out(y, x)*fscale) * (*W)(y, x);
                ans += delta * delta;
            }
        }
    }
    return ans;
}

template<class real>
real error_for_shift(Array<real> &out, const Array<double> &target, real &fscale) {
    fscale = filter_scale(out, target);
    real ans = error_given_fscale(out, target, fscale);
    if (is_symmetry()) {
        double t = profiling.current_iter * 1.0 / profiling.max_iterations;
        double symmetry_w = params.symmetry_w_min * pow(params.symmetry_w_max/params.symmetry_w_min, t);
        int w = target.width(), h = target.height();
        if (params.symmetry_h) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    real a = out(y, x);
                    real b = out(y, w-1-x);
                    real d = (a - b);
                    ans += symmetry_w * d * d;
                }
            }
        }
        if (params.symmetry_v) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    real a = out(y,     x);
                    real b = out(h-1-y, x);
                    real d = (a - b);
                    ans += symmetry_w * d * d;
                }
            }
        }
    }
    return ans;
}

template<class real>
real fit_error(const Array<real> &in, Array<real> &out, shared_ptr<Filter<real> > f, const Array<double> &target, bool rescale_out=false, real *fscale_out=NULL) {
    static vector<Array<real> > temp_images;
    static Array<real> in_p;
    static Array<real> out_p;
    static vector<Array<real> > out_L;
    
    real ans = 0.0, fscale = 0.0;
    if (!params.antialias || params.antialias_levels == 0) {
        f->apply({&in}, out, temp_images);
        ASSERT(out.sizes == target.sizes, "out size != K size");
        ans = error_for_shift(out, target, fscale);
        if (fscale_out) {
            *fscale_out = fscale;
        }
    } else {
        if (in.sizes != in_p.sizes) {
            in_p.resize(in.sizes);
        }
        int in_h = in.height(), in_w = in.width();
        vector<pair<int, int> > antialias_shift_L;
        int antialias_count = 0;
        
        /* End with (0, 0) so that the correct output data is written for the unshifted kernel */
        if (params.antialias_subsample <= 0) {
            antialias_shift_L.resize((1<<params.antialias_levels)*(1<<params.antialias_levels));
            for (int dy = (1<<params.antialias_levels)-1; dy >= 0; dy--) {
                for (int dx = (1<<params.antialias_levels)-1; dx >= 0; dx--) {
                    antialias_shift_L[antialias_count++] = pair<int, int>(dx, dy);
                }
            }
        } else {
            int level = MIN(params.antialias_levels, (int(antialias_subsample_x.size())-1));
            vector<int> &sub_x = antialias_subsample_x[level];
            vector<int> &sub_y = antialias_subsample_y[level];
            antialias_shift_L.resize(MIN(params.antialias_subsample, int(sub_x.size())));
            for (int j = 0; j < (int) antialias_shift_L.size(); j++) {
                antialias_shift_L[antialias_count++] = pair<int, int>(sub_x[j], sub_y[j]);
            }
        }
        out_L.resize(antialias_shift_L.size());
        for (int j = 0; j < (int) antialias_shift_L.size(); j++) {
            Array<real> &out_ref = j > 0 ? out_L[j]: out;
            int dx = antialias_shift_L[j].first;
            int dy = antialias_shift_L[j].second;
            for (int y = 0; y < in_h; y++) {
                int ysrc = y-dy;
                for (int x = 0; x < in_w; x++) {
                    int xsrc = x-dx;
                    if ((unsigned) xsrc < (unsigned) in_w && (unsigned) ysrc < (unsigned) in_h) {
                        in_p(y, x) = in(ysrc, xsrc);
                    } else {
                        in_p(y, x) = 0.0;
                    }
                }
            }
            f->apply({&in_p}, out_ref, temp_images);
            int out_h = out_ref.height(), out_w = out_ref.width();
            ASSERT(dx >= 0 && dy >= 0, "requires dx,dy>=0");
            for (int y = 0; y < out_h; y++) {
                int ysrc = y+dy;
                for (int x = 0; x < out_w; x++) {
                    int xsrc = x+dx;
                    if ((unsigned) xsrc < (unsigned) out_w && (unsigned) ysrc < (unsigned) out_h) {
                        out_ref(y, x) = out_ref(y+dy, x+dx);
                    } else {
                        out_ref(y, x) = 0.0;
                    }
                }
            }
            fscale += filter_scale(out_ref, target);
        }
        fscale /= antialias_count;
        for (int j = 0; j < (int) antialias_shift_L.size(); j++) {
            Array<real> &out_ref = j > 0 ? out_L[j]: out;
            ans += error_given_fscale(out_ref, target, fscale);
        }
        ans /= antialias_count;
        if (fscale_out) {
            *fscale_out = fscale;
        }
    }
    
    ASSERT(!params.search_shifts, "search_shifts unsupported");

    if (rescale_out) {
        for (int i = 0; i < target.nelems; i++) {
            out.data[i] *= fscale;
        }
    }
    real ans0 = ans;
#if PENALIZE_LOW_FSCALE
    real fscale_abs = fscale; //abs(fscale);
    if (fscale_abs < 0) { fscale_abs = -fscale_abs; }
    if (fscale_abs < 1) {
        real inv_fscale = 1.0/fscale_abs;
        ans += inv_fscale * LOW_FSCALE_PARAM;
    }
#else
    real fscale_abs = 0.0;
#endif
    if (params.verbose >= 3) {
        printf("fit_error: fscale=%e, err0(squared)=%f, fscale_abs=%f, err(squared)=%f\n", to_double(fscale), to_double(ans0), to_double(fscale_abs), to_double(ans));
        if (params.verbose >= 4) {
            for (int i = 0; i < target.nelems; i++) {
                printf("  target: %f, out: %f, out*fscale: %f\n", to_double(target.data[i]), to_double(out.data[i]), to_double(out.data[i]*fscale));
            }
            printf("\n");
        }
    }
    return ans;
}

template<class real>
double fit_error_double_sqrt(const Array<real> &in, Array<real> &out, shared_ptr<Filter<real> > f, const Array<double> &target, bool rescale_out=false) {
    RECLAIM_ADEPT_MEMORY();
    return qsqrt(to_double(fit_error(in, out, f, target, rescale_out)));
}

template<class real>
double fit_error_double_sqrt(const Array<real> &in, Array<real> &out, shared_ptr<FilterSparse<real> > f, const Array<double> &target, bool rescale_out=false) {
    return fit_error_double_sqrt(in, out, static_pointer_cast<Filter<real> >(f), target, rescale_out);
}

/* Numerical derivative (class double) */
/*
double fit_error_deriv(const Array<double> &in, Array<double> &out, shared_ptr<Filter<double> > f, const Array<double> &target, vector<double> &deriv);
*/

/* Automatic differentiation (class adouble) */
double fit_error_deriv(const Array<adouble> &in, Array<adouble> &out, shared_ptr<Filter<adouble> > f, const Array<double> &target, vector<double> &deriv);

template<class real>
Array<real> impulse(int w, int h) {
    vector<int> sizes({h, w});
    Array<real> ans(sizes);
    ans.clear();
    ans(h/2, w/2) = 1;
    return ans;
}

template<class real>
class BFGSParams { public:
    const Array<real> *in;
    Array<real> *out;
    shared_ptr<Filter<real> > f;
    const Array<double> *target;
    int niters;
    //bool verbose;
    vector<double> gradient;
    
    BFGSParams(const Array<real> *in_, Array<real> *out_, shared_ptr<Filter<real> > f_, const Array<double> *target_) :in(in_), out(out_), f(f_), target(target_), niters(0) {} //, verbose(verbose_) { }
};

template<class real>
static lbfgsfloatval_t lbfgs_evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    BFGSParams<real> *p = (BFGSParams<real> *) instance;
    vector<double> deriv;
    vector<real *> coeffs = p->f->coeffs();
    if (params.verbose >= 3) {
        printf("lbfgs_evaluate: coeffs:\n");
        for (int i = 0; i < n; i++) {
            printf("%f ", x[i]);
        }
        printf("\n");
    }
    for (int i = 0; i < n; i++) {
        (*coeffs[i]) = x[i];
    }
    double ans = fit_error_deriv(*p->in, *p->out, p->f, *p->target, deriv);

#if CHECK_NAN
    if (std::isnan(ans)) {
#if RECOVER_NAN
        for (int i = 0; i < (int) deriv.size(); i++) {
            g[i] = 0.0;
        }
        return 1e100;
#endif
        fprintf(stderr, "nan error\n");
        fprintf(stderr, "deriv: ");
        for (int i = 0; i < (int) deriv.size(); i++) {
            fprintf(stderr, "%f ", deriv[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "%s\n", p->f->str().c_str());
        fflush(stdout);
        fflush(stderr);
        assert(false);
        exit(1);
    }
#endif
    
    ASSERT(n == (int) deriv.size(), "n != deriv.size");
    //printf("deriv:\n%s\n\n", vector_to_str(deriv).c_str());
    //printf("%f\n", deriv[0]);
    
    if (params.verbose >= 2.5) {
        printf("lbfgs_evaluate: f=%f\n", ans);
        if (params.verbose >= 3) {
            printf("lbfgs_evaluate: gradient: \n");
            for (int i = 0; i < (int) deriv.size(); i++) {
                printf("%f ", deriv[i]);
            }
            printf("\n\n");
        }
    }
    for (int i = 0; i < (int) deriv.size(); i++) {
        g[i] = deriv[i];
#if CHECK_NAN
        if (std::isnan(deriv[i])) {
            fprintf(stderr, "nan derivative %d\n", i);
            fprintf(stderr, "%s\n", p->f->str().c_str());
            assert(false);
            exit(1);
        }
#endif
    }
    //printf("target: \n");
    //printf("%s\n\n", p->target->str().c_str());

    //printf("out: \n");
    //printf("%s\n\n", p->out->str().c_str());

    //printf("in: \n");
    //printf("%s\n\n", p->in->str().c_str());

    //printf("deriv: \n");
    //for (int i = 0; i < (int) deriv.size(); i++) {
    //    printf("%f ", deriv[i]);
    //}
    //printf("\n");
    
    //printf("ans=%f\n\n", ans);
    
    return ans;
}

template<class real>
static int lbfgs_progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
                    const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
    BFGSParams<real> *p = (BFGSParams<real> *) instance;
    p->niters = k;
    p->gradient.resize(n);
    for (int i = 0; i < n; i++) {
        p->gradient[i] = g[i];
    }
    if (params.verbose >= 3) {
        printf("Iteration %d: f = %f\n", k, fx);
        printf("Gradient: ");
        for (int i = 0; i < n; i++) {
            printf("%f ", p->gradient[i]);
        }
        printf("\n");
        printf("%s\n", p->f->str().c_str());
    }
    return 0;
}

template<class real>
double optimize_filter_bfgs(const Array<real> &in, Array<real> &out, shared_ptr<Filter<real> > f, const Array<double> &target0, int *iters=NULL, vector<double> *gradient=NULL, int solve_maxiters=DEFAULT_SOLVER_MAXITERS, double solve_epsilon=DEFAULT_SOLVER_EPSILON) {
    Array<double> target(target0);
    
    //printf("e\n"); fflush(stdout);
    //vector<double> deriv;
    //printf("e1\n"); fflush(stdout);
    BFGSParams<real> p(&in, &out, f, &target);
    //printf("e2\n"); fflush(stdout);
    auto coeffs = f->coeffs();
    //printf("e3\n"); fflush(stdout);
    int n = coeffs.size();
    //printf("e4 %d\n", n); fflush(stdout);
    lbfgsfloatval_t *x = new double[n];
    for (int i = 0; i < n; i++) {
        x[i] = to_double(*coeffs[i]);
    }
    //printf("e5\n"); fflush(stdout);
    lbfgsfloatval_t fx;
    //printf("f\n"); fflush(stdout);
    
    lbfgs_parameter_t lbfgs_params;
    lbfgs_parameter_init(&lbfgs_params);
    lbfgs_params.epsilon = solve_epsilon; //params.lbfgs_epsilon;
    lbfgs_params.max_iterations = solve_maxiters;
    int status = lbfgs(n, x, &fx, lbfgs_evaluate<real>, lbfgs_progress<real>, (void *) &p, &lbfgs_params);
    //printf("g\n"); fflush(stdout);

    for (int i = 0; i < n; i++) {
        *coeffs[i] = x[i];
    }

    if (params.verbose >= 2) {
        printf("L-BFGS optimization terminated with status=%d, %d iters (already minimized=%d)\n", status, p.niters, LBFGS_ALREADY_MINIMIZED);
        printf("  f = %f\n", fx);
        printf("\n");
        printf("%s\n", f->str().c_str());
    }
    
    if (iters) {
        *iters = p.niters;
    }
    if (gradient) {
        if (p.gradient.size() == 0) {
            lbfgsfloatval_t *g = new double[n];
            lbfgsfloatval_t step = 0, xnorm = 0, gnorm = 0;
            int k = 0, ls = 0;
            fx = lbfgs_evaluate<real>((void *) &p, x, g, n, step);
            lbfgs_progress<real>((void *) &p, x, g, fx, xnorm, gnorm, step, n, k, ls);
            delete[] g;
        }
        *gradient = p.gradient;
    }
    delete[] x;
    return fx;
}

template<class T>
T sum(const vector<T> &v) {
    T ans = 0;
    for (int i = 0; i < (int) v.size(); i++) {
        ans += v[i];
    }
    return ans;
}

#define RUN_BFGS(maxiters, epsilon) \
T0_bfgs = wall_time(); \
err = optimize_filter_bfgs(in, out, f_base, target, &bfgs_iters_read, &gradient, maxiters, epsilon); \
T1_bfgs = wall_time(); \
profiling.bfgs_calls++; \
profiling.bfgs_iters += bfgs_iters_read; \
profiling.T_bfgs += T1_bfgs - T0_bfgs;

/* Count zero coeffs in row/col, which has direction (span_dx, span_dy). Proceed to next row/col in direction (slab_dx, slab_dy). */
int count_zero_row_col(shared_ptr<Array<int> > box, const vector<MaskType> &mask, int x0, int y0, int span_dx, int span_dy, int slab_dx, int slab_dy);

vector<int> solve_translate(const vector<int> &dx_minL,const vector<int> &dx_maxL, const boxes_t &boxes);

template<class real>
void translate_box(shared_ptr<Array<int> > box, shared_ptr<FilterSparse<real> > f, int dx, int dy) {
    shared_ptr<FilterSparse<real> > f0 = f->copy_sameclass();
    
    int box_h = box->height();
    int box_w = box->width();
    vector<real *> f0_coeffs = f0->filter->coeffs();
    vector<real *> f_coeffs = f->filter->coeffs();
    ASSERT(f_coeffs.size() == f->mask.size(), "f_coeffs size != f->mask.size");
    for (int y = 0; y < box_h; y++) {
        for (int x = 0; x < box_w; x++) {
            int xp = x - dx;
            int yp = y - dy;
            MaskType mask_val = MASK_ZERO;
            real coeff_val = 0;
            if ((unsigned) xp < (unsigned) box_w && (unsigned) yp < (unsigned) box_h) {
                int idx_p = (*box)(yp, xp);
                ASSERT((unsigned) idx_p < (unsigned) f_coeffs.size(), "idx_p out of bounds");
                mask_val = f0->mask[idx_p];
                coeff_val = *f0_coeffs[idx_p];
            }
            
            int idx = (*box)(y, x);
            ASSERT((unsigned) idx < (unsigned) f_coeffs.size(), "idx out of bounds");
            f->mask[idx] = mask_val;
            *f_coeffs[idx] = coeff_val;
        }
    }
}

template<class real>
void find_convolution_chains(shared_ptr<FilterSparse<real> > f, vector<shared_ptr<boxes_t> > &chains) {
    //printf("find_convolution_chains\n");
    chains.clear();
    ASSERT(typeid(*f->filter) == typeid(FilterDAG<real>), "find_convolution_chains expected FilterSparse(FilterDAG())");
    shared_ptr<FilterDAG<real> > dag = dynamic_pointer_cast<FilterDAG<real> >(f->filter);

    map<int, int> parent_if_unique;                       // Parent index if has unique parent, -1 if not unique
    for (int i = 0; i < (int) dag->L.size(); i++) {
        for (int j = 0; j < (int) dag->L[i].inputs.size(); j++) {
            int input = dag->L[i].inputs[j];
            if (input != INPUT_SOURCE_IMAGE) {
                if (parent_if_unique.count(input)) {
                    parent_if_unique[input] = -1;
                } else {
                    parent_if_unique[input] = i;
                }
            }
        }
    }
    
    set<int> added;
    for (int i = 0; i < (int) dag->L.size(); i++) {
        if (added.count(i)) { continue; }
        int icurrent = i;
        boxes_t chain;
        while (1) {
            //printf("i = %d, icurrent = %d\n", i, icurrent);
            added.insert(icurrent);
            ASSERT((unsigned) icurrent < dag->L.size(), "icurrent out of bounds");
            shared_ptr<Filter<real> > current(dag->L[icurrent].f);
            bool continue_chain = false;
            if (typeid(*current) == typeid(FilterFIR<real>) || typeid(*current) == typeid(FilterFIRH<real>) || typeid(*current) == typeid(FilterFIRV<real>)) {
                //printf("  adding to chain FIR\n");
                boxes_t boxL = dag->subfilter_coeffs_boxes(icurrent);
                ASSERT(boxL.size() == 1, "expected 1 box for FilterIIR");
                chain.push_back(boxL[0]);
                continue_chain = true;
            } else if (typeid(*current) == typeid(FilterIIR<real>) ||
                       typeid(*current) == typeid(FilterIIRPlusY<real>) ||
                       typeid(*current) == typeid(FilterIIRMinusX<real>) ||
                       typeid(*current) == typeid(FilterIIRMinusY<real>)) {
                //printf("  adding to chain IIR\n");
                boxes_t boxL = dag->subfilter_coeffs_boxes(icurrent);
                ASSERT(boxL.size() == 2, "expected 2 boxes for FilterIIR");
                chain.push_back(boxL[0]);
                continue_chain = true;
            }
            //printf("continue_chain=%d\n", int(continue_chain));
            if (continue_chain) {
                if (parent_if_unique.count(icurrent)) {
                    //printf("parent_if_unique[icurrent=%d] = %d\n", icurrent, parent_if_unique[icurrent]);
                    if (parent_if_unique[icurrent] < 0) {
                        break;
                    } else {
                        ASSERT(parent_if_unique[icurrent] > icurrent, "parent invalid");
                        icurrent = parent_if_unique[icurrent];
                    }
                } else {
                    //printf("no parent_if_unique[icurrent=%d]\n", icurrent);
                    break;
                }
            } else {
                break;
            }
        }
        //printf("finished building chain of size %d\n", (int) chain.size());
        if (chain.size() >= 2) {
            chains.push_back(make_shared<boxes_t>(chain));
        }
    }
    //printf("done find_convolution_chains\n");
}

template<class real>
double translate_border(shared_ptr<FilterSparse<real> > f, const Array<real> &in, Array<real> &out, const Array<double> &target) {
    double err0 = 0.0;
#if DEBUG_TRANSLATE_BORDER
    shared_ptr<FilterSparse<real> > f0;
    boxes_t boxes_all;
    vector<int> dxL_all, dyL_all;
    vector<int> dx_maxL_all, dx_minL_all, dy_maxL_all, dy_minL_all;
#endif
    if (DEBUG_TRANSLATE_BORDER) {
        err0 = fit_error_double_sqrt<real>(in, out, dynamic_pointer_cast<Filter<real> >(f), target);
#if DEBUG_TRANSLATE_BORDER
        f0 = f->copy_sameclass();
#endif
    }
    
    vector<shared_ptr<boxes_t> > chains;
    find_convolution_chains(f, chains);
    
    for (int jchain = 0; jchain < (int) chains.size(); jchain++) {
        boxes_t boxes = *chains[jchain];
        
        if (boxes.size() != 2) { continue; }           // FIXME: Generalize to non 2 chains
        //ASSERT(typeid(dag->L[0].f) == typeid(FilterFIR<real>), "expected DAG node 0: FilterFIR");
        //ASSERT(typeid(dag->L[1].f) == typeid(FilterFIR<real>), "expected DAG node 1: FilterFIR");
        //shared_ptr<FilterFIR<real> > f0(dynamic_pointer_cast<FilterFIR<real> >(dag->L[0].f));
        //shared_ptr<FilterFIR<real> > f1(dynamic_pointer_cast<FilterFIR<real> >(dag->L[1].f));
        
        vector<int> dx_minL, dx_maxL;
        vector<int> dy_minL, dy_maxL;
        for (int i = 0; i < (int) boxes.size(); i++) {
            shared_ptr<Array<int> > box(boxes[i]);
            int w = box->width(), h = box->height();
            
            int zero_cols_left   = count_zero_row_col(box, f->mask, 0,   0, 0, 1,  1,  0);
            int zero_cols_right  = count_zero_row_col(box, f->mask, w-1, 0, 0, 1, -1,  0);
            int zero_rows_top    = count_zero_row_col(box, f->mask, 0,   0, 1, 0,  0,  1);
            int zero_rows_bottom = count_zero_row_col(box, f->mask, 0, h-1, 1, 0,  0, -1);
            
            dx_maxL.push_back(zero_cols_right-1);
            dx_minL.push_back(1-zero_cols_left);
            dy_maxL.push_back(zero_rows_bottom-1);
            dy_minL.push_back(1-zero_rows_top);
#if DEBUG_TRANSLATE_BORDER
            dx_maxL_all.push_back(zero_cols_right-1);
            dx_minL_all.push_back(1-zero_cols_left);
            dy_maxL_all.push_back(zero_rows_bottom-1);
            dy_minL_all.push_back(1-zero_rows_top);
#endif
        }
        
        vector<int> dxL = solve_translate(dx_minL, dx_maxL, boxes);
        vector<int> dyL = solve_translate(dy_minL, dy_maxL, boxes);
        
        for (int i = 0; i < (int) boxes.size(); i++) {
            translate_box<real>(boxes[i], f, dxL[i], dyL[i]);
        }
        
#if DEBUG_TRANSLATE_BORDER
        if (DEBUG_TRANSLATE_BORDER) {
            boxes_all.insert(boxes_all.end(), boxes.begin(), boxes.end());
            dxL_all.insert(dxL_all.end(), dxL.begin(), dxL.end());
            dyL_all.insert(dyL_all.end(), dyL.begin(), dyL.end());
        }
#endif
    }
    
    bool issue = false;
#if DEBUG_TRANSLATE_BORDER
    find_convolution_chains(f, chains);
    for (int jchain = 0; jchain < (int) chains.size(); jchain++) {
        boxes_t boxes = *chains[jchain];
        if (boxes.size() != 2) { continue; }           // FIXME: Generalize to non 2 chains
        
        for (int i = 0; i < (int) boxes.size(); i++) {
            shared_ptr<Array<int> > box(boxes[i]);
            int w = box->width(), h = box->height();
            
            int zero_cols_left   = count_zero_row_col(box, f->mask, 0,   0, 0, 1,  1,  0);
            int zero_cols_right  = count_zero_row_col(box, f->mask, w-1, 0, 0, 1, -1,  0);
            int zero_rows_top    = count_zero_row_col(box, f->mask, 0,   0, 1, 0,  0,  1);
            int zero_rows_bottom = count_zero_row_col(box, f->mask, 0, h-1, 1, 0,  0, -1);
            if (zero_cols_left == 0 || zero_cols_right == 0 || zero_rows_top == 0 || zero_rows_bottom == 0) {
                issue = true;
                break;
            }
        }
        if (issue) { break; }
    }
#endif

    double err1 = fit_error_double_sqrt<real>(in, out, dynamic_pointer_cast<Filter<real> >(f), target);
    if (DEBUG_TRANSLATE_BORDER) {
        double f_time = f->time();
        if (fabs(err1 - err0) > 1e-6) { // || (issue && f_time < 3.5 && err1 < 2)) {
            fprintf(stderr, "error changed after translate_border: %f, %f, or issue: %d\n", err0, err1, (int) issue);
            printf("time: %f\n", f_time);
            printf("mask:\n");
            for (int i = 0; i < (int) f->mask.size(); i++) {
                printf("%d ", f->mask[i]);
            }
            printf("\n");

            for (int jchain = 0; jchain < (int) chains.size(); jchain++) {
                boxes_t boxes = *chains[jchain];
                
                for (int i = 0; i < (int) boxes.size(); i++) {
                    shared_ptr<Array<int> > box(boxes[i]);
                    printf("chain %d, box %d\n", jchain, i);
                    for (int y = 0; y < box->height(); y++) {
                        for (int x = 0; x < box->width(); x++) {
                            printf("%d ", (*box)(y, x));
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }

#if DEBUG_TRANSLATE_BORDER
            for (int i = 0; i < (int) boxes_all.size(); i++) {
                fprintf(stderr, "box %d: %dx%d, dx,dy=%d,%d\n", i, boxes_all[i]->width(), boxes_all[i]->height(), dxL_all[i], dyL_all[i]);
            }
            for (int i = 0; i < (int) dx_maxL_all.size(); i++) {
                printf("dx/dy min_max: dx: %d %d dy: %d %d\n", dx_minL_all[i], dx_maxL_all[i], dy_minL_all[i], dy_maxL_all[i]);
            }
            printf("f0:\n%s\n\n", f0->str().c_str());
#endif
            printf("f:\n%s\n\n", f->str().c_str());
            ASSERT(false, "error changed after translate_border");
            exit(1);
        }
    }
    return err1;
}

template<class real>
void check_pareto(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto) {
    if (is_symmetry()) {
        return;
    }
    for (int i = 0; i < (int) pareto->L.size(); i++) {
        RECLAIM_ADEPT_MEMORY();
        double err = fit_error_double_sqrt(in, out, pareto->L[i].p, target);
        if (fabs(err - pareto->L[i].err) > 1e-6) {
            fprintf(stderr, "check_pareto: %d/%d, err=%f and pareto err=%f differ\n", i, (int) pareto->L.size(), err, pareto->L[i].err);
            assert(false);
        }
    }
}

string mask_to_str(const vector<MaskType> &m);

template<class real>
double optimize_and_add(const Array<real> &in, Array<real> &out, shared_ptr<FilterSparse<real> > f, const Array<double> &target, vector<double> &gradient, PARETO_TYPE(real) *pareto, bool *changed, int optimize_iter, shared_ptr<FilterSparse<real> > fparent) {
    if (params.mask_samples > 0) {
        string mask_str = mask_to_str(f->mask);
        if (profiling.mask_samples.count(mask_str) && profiling.mask_samples[mask_str] >= params.mask_samples) {
            profiling.skip_count++;
            if (profiling.skip_count % 10 == 0) { stack.new_recording(); }
            return OPTIMIZE_SKIP;
        }
        profiling.mask_samples[mask_str]++;
    }
    if (changed) { *changed = false; }
    int bfgs_iters_read = 0;
    //printf("mask before BFGS: %s\n", vector_to_str_int(f->mask).c_str());
    auto f_base = static_pointer_cast<Filter<real> >(f);
    double T0_bfgs, T1_bfgs, err;
    RUN_BFGS(params.solve_initial_maxiters, params.solve_initial_epsilon);
    if (params.verbose >= 1.9) {
        printf("Ran BFGS first time, got error: %f (before quantization)\n", sqrt(err));
    }
    profiling.num_solve_initial++;
    bool changed_v = false;
    if (!std::isnan(qsqrt(err))) {
        bool do_add = true;
        //printf("f after bfgs:\n%s\n\n", f->str().c_str());
        ////f->fix_mask();
        //printf("mask after fix_mask (1): %s\n", vector_to_str_int(f->mask).c_str());
        double f_time = f->time();
        double f_err = qsqrt(err);
        if (params.reoptim_add) {
            do_add = false;
            double dist = pareto->dist(f_time, f_err, f);
            if (dist <= params.reoptim_thresh) {
                RUN_BFGS(params.solve_add_maxiters, params.solve_add_epsilon);
                if (params.verbose >= 1.9) {
                    printf("Ran BFGS second time, got error: %f (before quantization)\n", sqrt(err));
                    fit_error_double_sqrt(in, out, f, target);
                }
                profiling.num_solve_add++;
                f_err = qsqrt(err);
                if (!std::isnan(f_err)) {
                    do_add = true;
                }
            }
        }
        //printf("do_add=%d\n", int(do_add));
        if (do_add) {
            f->fix_mask();
            if (params.translate_border) {
                //printf("mask after fix_mask (2): %s\n", vector_to_str_int(f->mask).c_str());
                f_err = translate_border(f, in, out, target);
            }
            auto fcopy = f->copy_sameclass();
            changed_v = pareto->add(f->time(), f_err, fcopy);
            if (params.adaptive) {
                profiling.adaptive_iter[(void *) (fcopy.get())] = optimize_iter;
                if (fparent) {
                    profiling.adaptive_iter[(void *) (fparent.get())] = optimize_iter;
                }
            }
            if (changed_v) {
                profiling.iter_changed = optimize_iter;
                profiling.add_count++;
                if (params.verbose >= 1.5) {
                    printf("Added Pareto time=%f, err=%f\n", f->time(), qsqrt(err));
                }
                if (params.check_often) {
                    check_pareto<real>(in, out, target, pareto);
                }
            }
            if (changed) {
                //check_pareto(in, out, target, pareto);
                *changed = changed_v;
            }
        }
    }
    if (params.verbose >= 2) {
        printf("  optimize_and_add: error %f, changed=%d\n", qsqrt(err), int(changed_v));
        printf("%s\n", f_base->str().c_str());
    }
    return err;
}

template<class real>
class FilterTuple { public:
    shared_ptr<FilterSparse<real> > f;
    double err;
    shared_ptr<FilterSparse<real> > fparent;
    FilterTuple(shared_ptr<FilterSparse<real> > f_, double err_, shared_ptr<FilterSparse<real> > fparent_) :f(f_), err(err_), fparent(fparent_) { }
};

template<class real>
bool operator < (const FilterTuple<real> &a, const FilterTuple<real> &b) {
    return a.err < b.err;
}

template<class real>
vector<int> get_num_samples(PARETO_TYPE(real) *pareto) {
    int n = (int) pareto->L.size();
    if (params.uniform_error) {
        vector<double> prob;
        //double Esum = 0.0;
        double Tmax = 0.0;
        for (int i = 0; i < n; i++) {
            Tmax = MAX(Tmax, pareto->L[i].T);
        }
        if (Tmax == 0.0) { Tmax = 1.0; }
        for (int i = 0; i < n; i++) {
            if (params.scale_correct && (pareto->L[i].err < 0 || pareto->L[i].err > 1+1e-6)) {
                fprintf(stderr, "pareto error < 0 or > 1: %f, index=%d\n", pareto->L[i].err, i);
                ASSERT(false, "err expected in [0, 1]");
                exit(1);
            }
            //ASSERT(pareto->L[i].err >= 0 && pareto->L[i].err <= 1, "err expected in [0, 1]");
            double Eprev = (i > 0) ? (fabs(pareto->L[i].err - pareto->L[i-1].err)): (2*MAX(1-pareto->L[i].err, 0));
            double Enext = (i+1 < n) ? (fabs(pareto->L[i+1].err - pareto->L[i].err)): (2*MAX(pareto->L[i].err, 0));
            double E = (Enext+Eprev)*0.5;

            double Tprev = (i > 0) ? (fabs(pareto->L[i].T/Tmax - pareto->L[i-1].T/Tmax)): (2*MAX(1-pareto->L[i].T/Tmax, 0));
            double Tnext = (i+1 < n) ? (fabs(pareto->L[i+1].T/Tmax - pareto->L[i].T/Tmax)): (2*MAX(pareto->L[i].T/Tmax, 0));
            double T = (Tnext+Tprev)*0.5;

            prob.push_back(sqrt(E*E+T*T));
            //Esum += E;
//            printf("err[%d] = %f, err[%d] = %f, err[%d] = %f, prob[%d] = %f, Esum = %f\n", i, pareto->L[i].err, i+1, (i+1 < n) ? pareto->L[i+1].err: -1, i-1, (i-1>0) ? pareto->L[i-1].err: -1, i, prob[i], Esum);
        }
//        printf("Esum = %f\n", Esum);
        ProbabilitySampler sampler(prob);
        
        vector<int> ans(n, 0);
        //printf("Building sample array: n=%d, len(prob)=%d\n", n, (int) prob.size());
        for (int i = 0; i < params.uniform_error_samples; i++) {
            int idx = sampler.sample();
        //    printf("  i=%d, idx=%d\n", i, idx);
            ans[idx]++;
        }
        
        return ans;
    } else {
        return vector<int>(n, 1);
    }
}

#define UPDATE_SAMPLES() \
    if (samples.size() != pareto->L.size()) { \
        samples = get_num_samples<real>(pareto); \
    } \

#define BEGIN_ITER_PARETO(start, cond, step) \
vector<int> samples; \
UPDATE_SAMPLES(); \
if (params.verbose >= 0.5) { \
    printf("Samples:\n"); \
    for (int ksample = 0; ksample < (int) samples.size(); ksample++) { \
        printf("%d ", samples[ksample]); \
    } \
    printf("\n"); \
} \
int jstep = (step); \
int jresample = 0; \
for (int j0 = (start); (cond); j0 += jstep) { \
    if (step < 0) { \
        j0 = MIN(j0, int(pareto->L.size())-1); \
        if (j0 < 0) { break; } \
    } \
    int j = params.random_order ? rand()%pareto->L.size(): j0; \
    jstep = (step); \
    if (params.adaptive) { \
        void *adaptive_ptr = (void *) (pareto->L[j].p.get()); \
        if (profiling.adaptive_iter.count(adaptive_ptr) && profiling.adaptive_iter[adaptive_ptr] < optimize_iter - params.adaptive_iters) { continue; } \
    } \
    UPDATE_SAMPLES(); \
    ASSERT((unsigned) j < samples.size(), "j out of bounds"); \
    for (int ksample = 0; ksample < samples[j]; ksample++) { \
        UPDATE_SAMPLES(); \
        if (j >= samples.size()) { break; }

#define CHECK_RESAMPLE(ferr) \
    if (ferr == OPTIMIZE_SKIP && jresample < params.max_resample) { jstep = 0; jresample++; continue; } \
    jresample = 0;

/*
    printf("j = %d, pareto size: %d, samples[%d] = %d\n", j, (int) pareto->L.size(), j, samples[j]); \
    for (int ii = 0; ii < samples.size(); ii++) { \
        printf("%d ", samples[ii]); \
    } \
    printf("\n"); \
        printf("j = %d, ksample = %d\n", j, ksample);
*/

#define END_ITER_PARETO() } }
    

/* Attempts to expand mask for each filter in the Pareto list. */
template<class real>
void optimize_sparse_filter_expand(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto, bool report_time=true, int optimize_iter=0) {
    profiling.reset();
    
    double T0 = wall_time();
    vector<double> gradient;
    
    BEGIN_ITER_PARETO(0, j0 < (int) pareto->L.size(), 1)
        if (params.verbose >= 1) {
            printf("optimize_sparse_filter_expand: j=%d/%d\n", j, (int) pareto->L.size());
        }
    
        //for (int k = 0; k < expand_iters; k++) {
        if (params.expand_is_random) {
            auto f = pareto->L[j].p->copy_sameclass();
            vector<int> mask_list;
            for (int i = 0; i < (int) f->mask.size(); i++) {
                if (f->mask[i] == MASK_ZERO) { mask_list.push_back(i); }
            }
            
            if (mask_list.size()) {
                int ichosen = mask_list[rand()%int(mask_list.size())];
                f->mask[ichosen] = MASK_VAR;
                //int orig_size = pareto->L.size();
                bool changed = false;
                double ferr = optimize_and_add(in, out, f, target, gradient, pareto, &changed, optimize_iter, pareto->L[j].p);
#if ATTEMPT_2MOVES
                if (ferr == OPTIMIZE_SKIP && jresample == params.max_resample) {
                    int ichosen2 = mask_list[rand()%int(mask_list.size())];
                    f->mask[ichosen2] = MASK_VAR;
                    ferr = optimize_and_add(in, out, f, target, gradient, pareto, &changed, optimize_iter, pareto->L[j].p);
                }
#endif
                CHECK_RESAMPLE(ferr);
                if (changed && params.expand_recheck) { jstep = 0; }
                //if (pareto->L.size() > orig_size) {
                //    j++;
                //}
            }
        } else {
            auto f0 = pareto->L[j].p;
            for (int i = 0; i < (int) f0->mask.size(); i++) {
//                    auto neighbors = mask_neighbors(i, f0->mask);
//                    for (int k = 0; k < (int) neighbors.size(); k++) {
//                        ASSERT(neighbors[k] >= 0 && neighbors[k] < (int) f0->mask.size(), "neighbor out of bounds");
//                        if (f0->mask[neighbors[k]] == MASK_ZERO) {
                if (f0->mask[i] == MASK_ZERO) {
                    auto f = f0->copy_sameclass();
                    f->mask[i] = MASK_VAR;
                    optimize_and_add(in, out, f, target, gradient, pareto, NULL, optimize_iter, f0);
//                        }
                }
            }
        }
        if (params.visit_one) { break; }
        //}
    END_ITER_PARETO();
    double T1 = wall_time();
    
    if (params.verbose || report_time) {
        printf("optimize_sparse_filter_expand (%d/%d): %f secs (%f BFGS, %d BFGS calls, %d BFGS iters, %f fit, %f deriv, %d OF calls, num_solve_initial: %d, num_solve_add: %d, added: %d, Pareto: %d)\n", optimize_iter, params.expand_iters, T1-T0, profiling.T_bfgs, profiling.bfgs_calls, profiling.bfgs_iters, profiling.T_fit, profiling.T_deriv, profiling.of_calls, profiling.num_solve_initial, profiling.num_solve_add, profiling.add_count, (int) pareto->L.size());
    }
}

template<class real>
void get_mask_bounding_box(shared_ptr<FilterSparse<real> > f, shared_ptr<Array<int> > box, int &row_min, int &row_max, int &col_min, int &col_max) {
    row_min = col_min = 1000*1000;
    row_max = col_max = 0;
    ASSERT(box->dimensions() == 2, "expected 2D box");
    for (int row = 0; row < box->height(); row++) {
        for (int col = 0; col < box->width(); col++) {
            int i = (*box)(row, col);
            ASSERT(i >= 0 && i < (int) f->mask.size(), "box(r, c) out of bounds");
            if (f->mask[i] == MASK_VAR) {
                row_min = MIN(row_min, row);
                row_max = MAX(row_max, row);
                col_min = MIN(col_min, col);
                col_max = MAX(col_max, col);
            }
        }
    }
}

template<class real>
void zero_row_col(vector<FilterTuple<real> > &L, shared_ptr<FilterSparse<real> > f, int row, int col, shared_ptr<Array<int> > box, const Array<real> &in, Array<real> &out, const Array<double> &target) {
    auto fcopy(f->copy_sameclass());
    for (int r = 0; r < box->height(); r++) {
        for (int c = 0; c < box->width(); c++) {
            if (row >= 0 && r != row) { continue; }
            if (col >= 0 && c != col) { continue; }
            int i = (*box)(r, c);
            ASSERT(i >= 0 && i < (int) fcopy->mask.size(), "box(r, c) out of bounds");
            fcopy->mask[i] = MASK_ZERO;
        }
    }
    double err_trunc = to_double(fit_error(in, out, static_pointer_cast<Filter<real> >(fcopy), target));
    L.push_back(FilterTuple<real>(fcopy, err_trunc, f));
    if (params.verbose >= 3) {
        printf("zero_row_col: %d %d success\n", row, col);
    }
}

int sample_box(boxes_t boxes);

/* Attempts to contract mask for each filter in the Pareto list. */
template<class real>
void optimize_sparse_filter_contract(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto, bool report_time=true, int optimize_iter=0) {
    double T0 = wall_time();
    double T_trunc = 0;
    double err = 1e100;
    vector<double> gradient;
    
    profiling.reset();
//    auto f(f0->copy_sameclass());
    if (params.verbose >= 1) {
        printf("optimize_sparse_filter_contract: begin\n");
    }

    BEGIN_ITER_PARETO(int(pareto->L.size()) - 1, j0 >= 0, -1)
        if (params.verbose >= 1) {
            printf("optimize_sparse_filter_contract: j=%d/%d\n", j, (int) pareto->L.size());
        }
        
        shared_ptr<FilterSparse<real> > f = pareto->L[j].p;
        if (params.truncate_taylor_approx) {
            err = fit_error_deriv(in, out, f, target, gradient);
        }
        
        int mask_count = 0;
        for (int i = 0; i < f->mask.size(); i++) {
            if (f->mask[i] != MASK_ZERO) { mask_count++; }
        }
        if (params.verbose >= 1) {
            printf("optimize_sparse_filter_contract: mask_count=%d\n", mask_count);
        }
/*        if (mask_count <= 1) {
            if (params.verbose >= 1) {
                printf("optimize_sparse_filter_contract: exiting due to mask_count=%d or gradient=%d\n", mask_count, (int) gradient.size());
            }
            break;
        }*/
        
        //double err_orig0 = to_double(fit_error(in, out, f.get(), target));
        double T0_trunc = wall_time();
        vector<FilterTuple<real> > L;
        auto coeffs = f->coeffs();
        //printf("coeffs size %d, gradient size %d\n", coeffs.size(), gradient.size());
        if (params.truncate_taylor_approx) {
            ASSERT(coeffs.size() == gradient.size(), "coeffs size != gradient size");
        }
        if (params.contract_is_random) {
            vector<int> sel_list;
            for (int i = 0; i < (int) f->mask.size(); i++) {
                if (f->mask[i] == MASK_VAR) {
                    sel_list.push_back(i);
                }
            }
            if (sel_list.size()) {
                int isel = sel_list[rand()%sel_list.size()];
                auto fcopy(f->copy_sameclass());
                fcopy->mask[isel] = MASK_ZERO;
                //auto fbase = static_pointer_cast<Filter<real> >(fcopy);
                bool changed = false;
                double ferr = optimize_and_add(in, out, fcopy, target, gradient, pareto, &changed, optimize_iter, f);
                CHECK_RESAMPLE(ferr);
#if ATTEMPT_2MOVES
                if (ferr == OPTIMIZE_SKIP && jresample == params.max_resample) {
                    int ichosen2 = sel_list[rand()%int(sel_list.size())];
                    f->mask[ichosen2] = MASK_ZERO;
                    ferr = optimize_and_add(in, out, f, target, gradient, pareto, &changed, optimize_iter, pareto->L[j].p);
                }
#endif
                if (changed && params.contract_recheck) {// && (j+1) < (int) pareto->L.size()) { //if (pareto->L.size() > old_size) {
                    jstep = 0;
                }
                continue;
                //L.push_back(FilterTuple<real>(fcopy, err_trunc, f));
            }
        } else {
            int ivar = 0;
            for (int i = 0; i < (int) f->mask.size(); i++) {
                //printf("mask variable %d/%d\n", i, (int) f->mask.size());
                if (f->mask[i] == MASK_VAR) {
                    auto fcopy(f->copy_sameclass());
                    fcopy->mask[i] = MASK_ZERO;
                    
                    double err_trunc = 0;
                    if (params.truncate_taylor_approx) {
                        ASSERT(ivar < (int) coeffs.size(), "ivar >= coeff size");
                        double dx = -to_double(*coeffs[ivar]);
                        err_trunc = err + dx * gradient[ivar];
                    } else {
                        //MaskType mask_orig = f->mask[i];
                        //f->mask[i] = MASK_ZERO;
                        auto fbase = static_pointer_cast<Filter<real> >(fcopy);
                        err_trunc = fit_error_double_sqrt<real>(in, out, fbase, target);  /* TODO: Could make this use double instead. Make a function for that. */
                        //f->mask[i] = mask_orig;
                    }
                    L.push_back(FilterTuple<real>(fcopy, err_trunc, f));
                    
                    ivar++;
                }
            }
        }
        if (params.try_separable) {
            boxes_t boxes(f->coeffs_boxes());
            int box_idx = sample_box(boxes);
            if (box_idx >= 0) {
                auto fcopy(f->copy_sameclass());
                int pattern = rand()%3;
                shared_ptr<Array<int> > box(boxes[box_idx]);
                
                int box_w = box->width(), box_h = box->height();
                int xvalid = -1, yvalid = -1;
                if (pattern == 0) {
                    xvalid = box_w/2;
                } else if (pattern == 1) {
                    yvalid = box_h/2;
                } else {
                    xvalid = box_w/2;
                    yvalid = box_h/2;
                }
                for (int y = 0; y < box_h; y++) {
                    for (int x = 0; x < box_w; x++) {
                        bool nonzero = true;
                        if (xvalid >= 0 && x != xvalid) { nonzero = false; }
                        if (yvalid >= 0 && y != yvalid) { nonzero = false; }
                        int isel = (*box)(y,x);
                        ASSERT((unsigned) isel < fcopy->mask.size(), "box index out of range");
                        fcopy->mask[isel] = nonzero ? MASK_VAR: MASK_ZERO;
                    }
                }
            //fcopy->mask[isel] = MASK_ZERO;
                auto fbase = static_pointer_cast<Filter<real> >(fcopy);
                L.push_back(FilterTuple<real>(fcopy, 0.0, f));
            }
        }
/*      // Buggy code -- this should make a copy of f instead.
        if (params.truncate_row_col && optimize_iter > 0) {
            boxes_t boxes(f->coeffs_boxes());
            //get_boxes_sparse(f, boxes);
            if (params.verbose >= 3) {
                printf("Number of boxes: %d\n", (int) boxes.size());
            }
            for (int ibox = 0; ibox < (int) boxes.size(); ibox++) {
                if (params.verbose >= 3) {
                    printf("ibox=%d\n", ibox);
                }
                shared_ptr<Array<int> > box(boxes[ibox]);
                int row_min = 0, row_max = 0, col_min = 0, col_max = 0;
                get_mask_bounding_box(f, box, row_min, row_max, col_min, col_max);
                if (params.verbose >= 3) {
                    printf("row_min=%d, row_max=%d, col_min=%d, col_max=%d\n", row_min, row_max, col_min, col_max);
                }
                if (row_max > row_min) {
                    zero_row_col(L, f, row_min, -1, box, in, out, target);
                    zero_row_col(L, f, row_max, -1, box, in, out, target);
                }
                if (col_max > col_min) {
                    zero_row_col(L, f, -1, col_min, box, in, out, target);
                    zero_row_col(L, f, -1, col_max, box, in, out, target);
                }
            }
        }
*/
        double T1_trunc = wall_time();
        T_trunc += T1_trunc - T0_trunc;
        
        sort(L.begin(), L.end());
        
        double err_best = 1e100;
        bool success = false;
        for (int i = 0; i < MIN(params.optimize_top_k, (int) L.size()); i++) {
            //int old_size = pareto->L.size();
            //auto old_L = vector<ParetoPoint<FilterSparse<real> > >(pareto->L);
            bool changed = false;
            double err_opt = optimize_and_add(in, out, L[i].f, target, gradient, pareto, &changed, optimize_iter, L[i].fparent);
            if (changed && params.contract_recheck) {// && (j+1) < (int) pareto->L.size()) { //if (pareto->L.size() > old_size) {
                jstep = 0;
            }
            if (err_opt < err_best) {
                err_best = err_opt;
                f = L[i].f;
                success = true;
            }
        }
        err = err_best;
        if (!success) {
            if (params.verbose >= 1) {
                printf("No success finding truncated filter, stopping\n");
            }
            break;
        }
        if (params.visit_one) { break; }
    END_ITER_PARETO();
    double T1 = wall_time();

    if (params.verbose || report_time) {
        printf("optimize_sparse_filter_contract (%d/%d): %f secs (%f BFGS, %d BFGS calls, %d BFGS iters, %f fit, %f deriv, %d OF calls, num_solve_initial: %d, num_solve_add: %d, added: %d, Pareto: %d, mask_samples: %d)\n", optimize_iter, params.contract_iters, T1-T0, profiling.T_bfgs, profiling.bfgs_calls, profiling.bfgs_iters, profiling.T_fit, profiling.T_deriv, profiling.of_calls, profiling.num_solve_initial, profiling.num_solve_add, profiling.add_count, (int) pareto->L.size(), (int) profiling.mask_samples.size());
    }
}

template<class real>
void optimize_sparse_filter_translate(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto, bool report_time=true, int optimize_iter=0) {
    profiling.reset();
    
    double T0 = wall_time();
    vector<double> gradient;
    
    int repeats = 1;
    if (params.translate_step < 1) { repeats = int(1.0/params.translate_step); }
    for (int repeat = 0; repeat < repeats; repeat++) {
        ASSERT(pareto->L.size(), "pareto list empty");
        for (int kiter = 0; kiter < 5; kiter++) {
            bool do_all = params.translate_all && kiter == 0;
            bool do_some = params.translate_some && kiter == 1;
            bool do_duplicate = params.translate_duplicate && kiter == 2;
            bool do_cross = params.translate_cross && kiter == 3 && (pareto->L.size() >= 2);
            bool do_copy = params.translate_copy && kiter == 4;
            if (!do_all && !do_some && !do_duplicate && !do_cross && !do_copy) { continue; }

            BEGIN_ITER_PARETO(0, j0 < (int) pareto->L.size(), 1)
                if (params.verbose >= 1) {
                    printf("optimize_sparse_filter_translate: j=%d/%d\n", j, (int) pareto->L.size());
                }

                auto f0 = pareto->L[j].p;
            
                auto f = f0->copy_sameclass();
                boxes_t boxes;
                int box_idx = -1;
                shared_ptr<Array<int> > box;
                
                if (do_all || do_some || do_duplicate || do_copy) {
                    boxes = f->coeffs_boxes();
                    box_idx = sample_box(boxes);
                    if (box_idx < 0) { continue; }
                    box = boxes[box_idx];
                }
            
                if (do_copy) {        /* Copy between boxes of identical size */
                    //printf("do_copy, box_idx=%d\n", box_idx);
                    vector<int> src_idx_L;
                    for (int i = 0; i < (int) boxes.size(); i++) {
                        if (boxes[i]->sizes == box->sizes && i != box_idx) {
                            src_idx_L.push_back(i);
                        }
                    }
                    if (src_idx_L.size()) {
                        int src_idx = src_idx_L[rand()%src_idx_L.size()];
                        //printf("copying from %d -> %d\n", box_idx, src_idx);
                        shared_ptr<Array<int> > src_box = boxes[src_idx];
                        
                        ASSERT(src_box->sizes == box->sizes, "src_box and box differ in size");

                        int box_h = box->height();
                        int box_w = box->width();
                        ASSERT(box_w > 0 && box_h > 0, "box dims should be > 0");
                        
                        for (int y = 0; y < box_h; y++) {
                            for (int x = 0; x < box_w; x++) {
                                f->mask[(*box)(y,x)] = f->mask[(*src_box)(y,x)];
                            }
                        }
                    }
                }
            
                if (do_all) {         /* Translate entire box */
                    int dir = rand()%4;
                    int dx = 0, dy = 0;
                    
                    if (dir == 0) { dx = -1; }
                    else if (dir == 1) { dx = 1; }
                    else if (dir == 2) { dy = -1; }
                    else if (dir == 3) { dy = 1; }
                    else { ASSERT(false, "invalid case"); }
                    
                    translate_box<real>(box, f, -dx, -dy);
                }

                if (do_some) {         /* Translate individual elements */
                    vector<int> adj_idx0, adj_idx1;
                    ASSERT(box->dimensions() == 2, "box should be 2D");
                    int box_h = box->height();
                    int box_w = box->width();
                    for (int y = -1; y < box_h; y++) {
                        for (int x = -1; x < box_w; x++) {
                            for (int di = 0; di < 2; di++) {
                                int dx = di == 0;
                                int dy = di == 1;
                                int xp = x + dx;
                                int yp = y + dy;
                                bool bounds0 = (unsigned) x  < (unsigned) box_w && (unsigned) y  < (unsigned) box_h;
                                bool bounds1 = (unsigned) xp < (unsigned) box_w && (unsigned) yp < (unsigned) box_h;
                                bool mask0 = false, mask1 = false;
                                int idx0 = -1, idx1 = -1;
                                if (bounds0) {
                                    idx0 = (*box)(y,x);
                                    ASSERT((unsigned) idx0 < f->mask.size(), "idx0 out of bounds");
                                    mask0 = f->mask[idx0] == MASK_VAR;
                                }
                                if (bounds1) {
                                    idx1 = (*box)(yp,xp);
                                    ASSERT((unsigned) idx1 < f->mask.size(), "idx1 out of bounds");
                                    mask1 = f->mask[idx1] == MASK_VAR;
                                }
                                if (mask0 ^ mask1) {
                                    adj_idx0.push_back(idx0);
                                    adj_idx1.push_back(idx1);
                                }
                            }
                        }
                    }

                    ASSERT(adj_idx0.size(), "adj_idx0 empty");
                    int iadj = rand()%adj_idx0.size();
                    int idx0 = adj_idx0[iadj];
                    int idx1 = adj_idx1[iadj];
                    ASSERT(idx0 >= 0 || idx1 >= 0, "both indices < 0");
                    bool mask0 = idx0 >= 0 ? f->mask[idx0] == MASK_VAR: false;
                    bool mask1 = idx1 >= 0 ? f->mask[idx1] == MASK_VAR: false;
                    
                    bool mask0p = false, mask1p = false;
                    while (1) {
                        mask0p = bool(rand()%2);
                        mask1p = bool(rand()%2);
                        if (idx0 < 0) { mask0p = mask0; }
                        if (idx1 < 0) { mask1p = mask1; }
                        if (mask0p != mask0 || mask1p != mask1) { break; }
                    }

                    if (idx0 >= 0) { f->mask[idx0] = mask0p ? MASK_VAR: MASK_ZERO; }
                    if (idx1 >= 0) { f->mask[idx1] = mask1p ? MASK_VAR: MASK_ZERO; }
                }
                
                if (do_duplicate) {         /* Duplicate a row or column */
                    ASSERT(box->dimensions() == 2, "box should be 2D");
                    int box_h = box->height();
                    int box_w = box->width();
                    ASSERT(box_w > 0 && box_h > 0, "box dims should be > 0");
                    
                    int dup_row = -1, dup_col = -1;
                    int dx = 0, dy = 0;
                    if (rand() % 2 == 0) {
                        dup_row = params.translate_duplicate_outside ? (-1+(rand()%(box_h+2))): (rand()%box_h);
                        dy = rand() % 2 == 0 ? -1: 1;
                    } else {
                        dup_col = params.translate_duplicate_outside ? (-1+(rand()%(box_w+2))): (rand()%box_w);
                        dx = rand() % 2 == 0 ? -1: 1;
                    }
                    
                    for (int y = 0; y < box_h; y++) {
                        for (int x = 0; x < box_w; x++) {
                            int xsrc = x, ysrc = y;
                            if (dx > 0) {
                                if (x > dup_col) { xsrc--; }
                            } else if (dx < 0) {
                                if (x < dup_col) { xsrc++; }
                            } else if (dy > 0) {
                                if (y > dup_row) { ysrc--; }
                            } else if (dy < 0) {
                                if (y < dup_col) { ysrc++; }
                            } else {
                                ASSERT(false, "invalid dx, dy");
                            }
                            int idx = (*box)(y,x);
                            if ((unsigned) xsrc < (unsigned) box_w && (unsigned) ysrc < (unsigned) box_h) {//, "xsrc,ysrc out of bounds");
                                int idx_src = (*box)(ysrc, xsrc);
                                ASSERT((unsigned) idx_src < f->mask.size(), "idx_src out of bounds");
                                ASSERT((unsigned) idx < f->mask.size(), "idx out of bounds");
                                f->mask[idx] = f->mask[idx_src];
                            } else {
                                f->mask[idx] = MASK_ZERO;
                            }
                        }
                    }
                }

                if (do_cross) {
                    int idx2 = j;
                    while (idx2 == j) {
                        idx2 = rand()%pareto->L.size();
                    }
                    shared_ptr<FilterSparse<real> > f2(pareto->L[idx2].p);
                    vector<real *> f_coeffs(f->filter->coeffs());
                    vector<real *> f2_coeffs(f2->filter->coeffs());
                    int cut = 1;
                    if (f_coeffs.size() > 2) {
                        cut = 1+rand()%(f_coeffs.size()-1);
                    }
                    
                    for (int i = cut; i < (int) f2_coeffs.size(); i++) {
                        *f_coeffs[i] = *f2_coeffs[i];
                    }
                    f->mask_from_coeffs();
                }
                
                CHECK_RESAMPLE(optimize_and_add(in, out, f, target, gradient, pareto, NULL, optimize_iter, f));
                if (params.visit_one) { break; }
            END_ITER_PARETO();
        }
    }
    
    double T1 = wall_time();
    
    if (params.verbose || report_time) {
        printf("optimize_sparse_filter_translate: %f secs (%f BFGS, %d BFGS calls, %d BFGS iters, %f fit, %f deriv, %d OF calls)\n", T1-T0, profiling.T_bfgs, profiling.bfgs_calls, profiling.bfgs_iters, profiling.T_fit, profiling.T_deriv, profiling.of_calls);
    }
}

void detect_symmetry(const Array<double> &target);

template<class real>
void update_pareto(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto) {
    PARETO_TYPE(real) pareto_copy(*pareto);
    pareto->clear();
    for (int i = 0; i < (int) pareto_copy.L.size(); i++) {
        shared_ptr<Filter<real> > fbase = dynamic_pointer_cast<Filter<real> >(pareto_copy.L[i].p);
        double err = fit_error_double_sqrt(in, out, fbase, target);
        pareto->add(pareto_copy.L[i].T, err, pareto_copy.L[i].p);
    }
}

template<class real>
void randomize_filter(shared_ptr<FilterSparse<real> > fcopy, int num_nonzero=-1, bool restrict_values=false) {
    int ncoeffs = fcopy->coeffs().size();
    if (num_nonzero < 0) {
        num_nonzero = rand()%(ncoeffs+1);
    }
    vector<real *> coeffs = fcopy->filter->coeffs();
    vector<int> idx_nonzero;
    for (int j = 0; j < ncoeffs; j++) {
        idx_nonzero.push_back(j);
    }
    std::random_shuffle(idx_nonzero.begin(), idx_nonzero.end());
    int allow_negative = rand()%4 == 0;
    for (int j = 0; j < ncoeffs; j++) {
        fcopy->mask[j] = MASK_ZERO;
        *coeffs[j] = 0.0;
    }
    for (int j = 0; j < num_nonzero; j++) {
        int idx = idx_nonzero[j];
        fcopy->mask[idx] = MASK_VAR;
        if (restrict_values) {
            double maxv = 1.0/100.0;
            *coeffs[idx] = rand_f() * maxv;
        } else {
            *coeffs[idx] = allow_negative ? (2*rand_f()-1): rand_f();
        }
    }
    fcopy->mask_from_coeffs();
}

template<class real>
void inject_random(const Array<real> &in, Array<real> &out, shared_ptr<FilterSparse<real> > f0, const Array<double> &target, PARETO_TYPE(real) *pareto) {
    if (!params.random_init) { return; }
    vector<double> gradient;
    int ncoeffs = f0->coeffs().size();
    int mask_samples0 = params.mask_samples;
    int random_count0 = params.random_count;
    int solve_initial_maxiters0 = params.solve_initial_maxiters;
    params.mask_samples = -1;
    //if (params.antialias) {
        //if (params.antialias_levels >= 2) {
        //    params.random_count = MIN(params.random_count, 50);
    //    params.solve_initial_maxiters = 100;
        //}
    //}
    for (int i = 0; i < params.random_count; i++) {
        int num_nonzero = ncoeffs >= 1 ? (rand()%(ncoeffs+1)): ncoeffs;  // FIXME: Should we disable mask_samples for inject_random?
        if (i == 0) {
            num_nonzero = ncoeffs;
            params.solve_initial_maxiters = 2000;
        }
        if (params.verbose >= 1.9) { printf("inject_random: %d/%d, num_nonzero=%d\n", i, params.random_count, num_nonzero); }
        
        shared_ptr<FilterSparse<real> > fcopy(f0->copy_sameclass());
        randomize_filter(fcopy, num_nonzero);
#if DEBUG_MULTIBAND
        boxes_t boxes = fcopy->coeffs_boxes();
        printf("DEBUG_MULTIBAND, boxes: %d\n", (int) boxes.size());
        if (boxes.size() == 7) {
            int mask_box_min = 3, mask_box_max = 5;
            for (int mask_box = mask_box_min; mask_box <= mask_box_max; mask_box++) {
                printf("DEBUG_MULTIBAND, %d boxes, masking out box %d", (int) boxes.size(), mask_box);
                auto box = boxes[mask_box];
                for (int i = 0; i < box->nelems; i++) {
                    int idx = (*box).data[i];
                    ASSERT2((unsigned) idx < (unsigned) coeffs.size(), "idx out of bounds");
                    *coeffs[idx] = 0.0;
                    fcopy->mask[idx] = MASK_ZERO;
                }
            }
        }
#endif
        optimize_and_add(in, out, fcopy, target, gradient, pareto, NULL, 0, SPARSE_NULL);
        if (i == 0) {
            params.solve_initial_maxiters = solve_initial_maxiters0;
        }
    }
    params.mask_samples = mask_samples0;
    params.random_count = random_count0;
    params.solve_initial_maxiters = solve_initial_maxiters0;
}

/* Get solution for just one FIR filter */
template<class real>
void get_base_solution(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto) {
//    vector<int> K_sizes({params.filter_w, params.filter_w});
    Array<real> K(target.sizes);
    for (int y = 0; y < target.height(); y++) {
        for (int x = 0; x < target.width(); x++) {
//            int xp = (x-params.filter_w/2) + target.width()/2;
//            int yp = (y-params.filter_w/2) + target.height()/2;
//            if ((unsigned) xp < (unsigned) target.width() && (unsigned) yp < (unsigned) target.height()) {
//                K(y, x) = target(yp, xp);
//            } else {
//                K(y, x) = 0.0;
//            }
            K(y, x) = target(y, x);
        }
    }
    
    auto f0_fir = make_shared<FilterFIR<real> >(K);
    auto f0_fir_base = static_pointer_cast<Filter<real> >(f0_fir);
    auto f0_fir_baseL = vector<shared_ptr<Filter<real> > >({f0_fir_base});
    
    auto f0_dag = make_shared<FilterDAG<real> >(f0_fir_baseL);
    auto f0 = make_shared<FilterSparse<real> >(static_pointer_cast<Filter<real> >(f0_dag));
    f0->mask_from_coeffs();
    auto f0_base = static_pointer_cast<Filter<real> >(f0);

    double f_err = fit_error_double_sqrt(in, out, f0_base, target);
    pareto->add(f0->time(), f_err, f0->copy_sameclass());
    while (1) {
        vector<real *> f_coeffs = f0->filter->coeffs();
        ASSERT(f_coeffs.size() == f0->mask.size(), "coeffs size differs from mask size");
        double coeff_least_abs = 1e100;
        int coeff_least_idx = -1;
        for (int i = 0; i < (int) f0->mask.size(); i++) {
            if (f0->mask[i] == MASK_ZERO) { continue; }
            double c = to_double(*f_coeffs[i]);
            if (c < 0) { c = -c; }
            if (c < coeff_least_abs) {
                coeff_least_abs = c;
                coeff_least_idx = i;
            }
        }
        if (coeff_least_idx < 0) { break; }
        f0->mask[coeff_least_idx] = MASK_ZERO;
        *f_coeffs[coeff_least_idx] = 0.0;
        
        f_err = fit_error_double_sqrt(in, out, f0_base, target);
        pareto->add(f0->time(), f_err, f0->copy_sameclass());
    }
}

template<class real>
Array<real> multisize_resize_array(const Array<real> &K, int sz, const Array<real> &K0) {
    int w = K.width();
    int h = K.height();
    int w0 = K0.width();
    int h0 = K0.height();
    int wp = MIN(sz, w0);
    int hp = MIN(sz, h0);
    
    vector<int> sizes({hp, wp});
    Array<real> ans(sizes);
    for (int y = 0; y < hp; y++) {
        for (int x = 0; x < wp; x++) {
            int xsrc = x-wp/2+w/2;
            int ysrc = y-hp/2+h/2;
            if ((unsigned) xsrc < (unsigned) w && (unsigned) ysrc < (unsigned) h) {
                ans(y, x) = K(ysrc, xsrc);
            } else {
                ans(y, x) = 0.0;
            }
        }
    }
    return ans;
}

template<class real>
shared_ptr<FilterSparse<real> > multisize_resize_filter(shared_ptr<FilterSparse<real> > f, int sz, shared_ptr<FilterSparse<real> > f0) {
    f = f->copy_sameclass();
    f->mask_from_coeffs();
    ASSERT(typeid(*f->filter) == typeid(FilterDAG<real>), "expected FilterSparse(FilterDAG())");
    shared_ptr<FilterDAG<real> > f_dag(static_pointer_cast<FilterDAG<real> >(f->filter));
    shared_ptr<FilterDAG<real> > f0_dag(static_pointer_cast<FilterDAG<real> >(f0->filter));

//    boxes_t boxes0(f0->coeffs_boxes());
//    boxes_t boxes(f->coeffs_boxes());
//    ASSERT(boxes0.size() == boxes.size(), "expected f0 and f to have same boxes size");
//    ASSERT(f_dag->L.size() == boxes.size(), "expected DAG size to be boxes size");
    for (int i = 0; i < (int) f_dag->L.size(); i++) {
        if (typeid(*f_dag->L[i].f) == typeid(FilterFIR<real>)) {
            auto fcurrent  = static_pointer_cast<FilterFIR<real> >(f_dag->L[i].f);
            auto fcurrent0 = static_pointer_cast<FilterFIR<real> >(f0_dag->L[i].f);
            fcurrent->K = multisize_resize_array(fcurrent->K, sz, fcurrent0->K);
        } else if (typeid(*f_dag->L[i].f) == typeid(FilterIIR<real>)) {
            auto fcurrent  = static_pointer_cast<FilterIIR<real> >(f_dag->L[i].f);
            auto fcurrent0 = static_pointer_cast<FilterIIR<real> >(f0_dag->L[i].f);
            fcurrent->K = multisize_resize_array(fcurrent->K, sz, fcurrent0->K);
            fcurrent->F = multisize_resize_array(fcurrent->F, sz, fcurrent0->F);
        }
    }
    f->mask_from_coeffs();
    return f;
}

template<class real>
void multisize_resize_pareto(PARETO_TYPE(real) *pareto, int sz, shared_ptr<FilterSparse<real> > f0, const Array<real> &in, Array<real> &out, const Array<double> &target) {
    PARETO_TYPE(real) pareto0(*pareto);
    pareto->clear();
    for (int i = 0; i < (int) pareto0.L.size(); i++) {
        shared_ptr<FilterSparse<real> > f(pareto0.L[i].p);
        double T0 = pareto0.L[i].T;
        double err0 = pareto0.L[i].err;
        f = multisize_resize_filter(f, sz, f0);
        double T = f->time();
        double err = fit_error_double_sqrt(in, out, f, target);
        ASSERT(fabs(T0-T) < 1e-7, "T and T0 differ");
        ASSERT(fabs(err0-err) < 1e-7, "err and err0 differ");
        pareto->add(T, err, f);
    }
}

vector<int> get_multisize_sizes();

bool terminate_early(const vector<double> &Tavg_L, double T_start_optim);

#define STOP_IF_DONE() \
if (params.max_of_calls > 0 && profiling.total_of_calls - start_of_calls > params.max_of_calls/multisize_sizes.size()) { break; }

template<class real>
int dag_levels(shared_ptr<FilterDAG<real> > dag, vector<double> &sizes) {
    if (dag->L.size() == 0) { return 0; }
    dag->check_sizes(sizes);
    if (sizes.size() == 0) { return 0; }
    
    double min_size = *std::min_element(sizes.begin(), sizes.end());
    int levels = 0;
    while (min_size < 1-1e-8) {
        min_size *= 4;
        levels++;
    }
    return levels;
}

//template<class real>
//void set_antialias_levels(shared_ptr<FilterSparse<real> > f0) {
//    ASSERT2(typeid(*f0->filter) == typeid(FilterDAG<real>), "expected FilterSparse(FilterDAG()) in set_antialias_levels");
//    shared_ptr<FilterDAG<real> > dag = dynamic_pointer_cast<FilterDAG<real> >(f0->filter);
//    vector<double> sizes;
//    params.antialias_levels = dag_levels(dag, sizes);
//    printf("set_antialias_levels %d\n", params.antialias_levels);
//}

/* Iteratively sparsifies filter and optimizes remaining coeffs using BFGS, adding each result to the Pareto frontier */
template<class real>
void optimize_sparse_filter(const Array<real> &in, Array<real> &out, shared_ptr<FilterSparse<real> > f0, const Array<double> &target, PARETO_TYPE(real) *pareto, bool report_time=true, bool init_profiling=true) {
    bool quantize0 = params.quantize;
    double T_start_optim = wall_time();
    //if (params.antialias) {
    //    set_antialias_levels(f0);
    //}
    if (init_profiling) {
        profiling.init();
    }
    vector<double> gradient;
    vector<double> Tavg_L;

    /* If no tunable parameters, return single solution. */
    if (f0->coeffs().size() == 0) {
        auto fcopy = f0->copy_sameclass();
        double f_err = fit_error_double_sqrt(in, out, fcopy, target);
        pareto->add(fcopy->time(), f_err, fcopy);
        write_pareto_trace(pareto, params.trace, T_start_optim, profiling.total_of_calls, in, out, target, &Tavg_L);
        return;
    }
    
    vector<int> multisize_sizes = get_multisize_sizes();
    
    params.allow_symmetry = true;
    detect_symmetry(target);
    
    bool terminated = false;
    int iter = 0;
    for (int imultisize = 0; imultisize < (int) multisize_sizes.size(); imultisize++) {
        if (params.quantize_until >= 0 && params.quantize && iter >= params.quantize_until) {
            params.quantize = false;
            update_pareto(in, out, target, pareto);
        }
        int sz = multisize_sizes[imultisize];
        int start_of_calls = profiling.total_of_calls;
        
        shared_ptr<FilterSparse<real> > f0resized(f0);
        if (params.multisize) {
            f0resized = multisize_resize_filter(f0, sz, f0);
            multisize_resize_pareto(pareto, sz, f0, in, out, target);
            //printf("----------------- Pareto -----------------\n");
            //for (int i = 0; i < (int) pareto->L.size(); i++) {
            //    printf("%s\n\n", pareto->L[i].p->str().c_str());
            //}
            //printf("------------------------------------------\n");
        }
        if (params.verbose || report_time) {
            printf("=== Multisize mode=%d imultisize=%d/%d, current size=%d ===\n", params.multisize, imultisize, (int) multisize_sizes.size(), sz);
        }
        
        optimize_and_add(in, out, f0resized->copy_sameclass(), target, gradient, pareto, NULL, 0, SPARSE_NULL);
        
        if (params.resume.size() == 0) {
            inject_random(in, out, f0resized->copy_sameclass(), target, pareto);
        }
        
        int iters = params.expand_iters;
        iters = MAX(iters, params.contract_iters);
        iters = MAX(iters, params.translate_iters);
        profiling.iter_changed = 0;
        profiling.max_iterations = iters;
        
        int last_iter = iter+iters/multisize_sizes.size();
        for (; iter < last_iter; iter++) {
            profiling.current_iter = iter;
            if (is_symmetry()) {
                if (iter % params.symmetry_reset_iters == 0 && params.mask_samples > 0) {
                    profiling.mask_samples.clear();
                }
                update_pareto(in, out, target, pareto);
            }
            if (iter < params.contract_iters) {
                optimize_sparse_filter_contract(in, out, target, pareto, report_time, iter);
            }
            STOP_IF_DONE();
            //check_pareto(in, out, target, pareto);
            if (iter < params.expand_iters) {
                optimize_sparse_filter_expand(in, out, target, pareto, report_time, iter);
            }
            STOP_IF_DONE();
            //check_pareto(in, out, target, pareto);
            int translate_step_ipart = MAX(int(params.translate_step), 1);
            if (iter < params.translate_iters && iter % translate_step_ipart == 0) {
                optimize_sparse_filter_translate(in, out, target, pareto, report_time, iter);
            }
            STOP_IF_DONE();
            //check_pareto(in, out, target, pareto);
            if (params.verbose || report_time) {
                fflush(stdout);
            }
            if (is_symmetry()) {
                params.allow_symmetry = false;
                update_pareto(in, out, target, pareto);
                params.allow_symmetry = true;
            }
            write_pareto_trace(pareto, params.trace, T_start_optim, profiling.total_of_calls, in, out, target, &Tavg_L);
            if (terminate_early(Tavg_L, T_start_optim)) {
                terminated = true;
                break;
            }
        }
        if (terminated) { break; }
    }
    
    if (params.verbose || report_time) {
        printf("Last changed on iteration %d\n", profiling.iter_changed);
    }
    params.quantize = quantize0;
    
//    if (params.trace.size()) {
//        write_pareto_trace(pareto, params.trace, true);
//    }
}

#endif

