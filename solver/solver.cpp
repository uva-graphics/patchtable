
#include "solver.h"

adept::Stack stack;

/* Antialiasing (dx, dy) patterns output by 'python hammingseq.py 9'. */
vector<vector<int> > antialias_subsample_x({{0}, {0, 1, 0, 1}, {0, 3, 0, 1, 1, 2, 2, 3}, {0, 7, 0, 3, 5, 6, 1, 1}, {0, 15, 0, 3, 3, 5, 5, 6}, {0, 31, 0, 7, 11, 21, 26, 28}, {0, 63, 0, 7, 7, 25, 25, 30}, {0, 127, 0, 15, 23, 57, 90, 99}, {0, 255, 0, 15, 15, 51, 51, 60}, {0, 511, 0, 31, 47, 211, 220, 355}});
vector<vector<int> > antialias_subsample_y({{0}, {0, 1, 1, 0}, {0, 3, 3, 1, 2, 1, 2, 0}, {0, 7, 7, 1, 2, 4, 3, 4}, {0, 15, 15, 3, 12, 5, 10, 6}, {0, 31, 31, 3, 12, 20, 17, 10}, {0, 63, 63, 7, 56, 11, 52, 12}, {0, 127, 127, 7, 56, 73, 82, 100}, {0, 255, 255, 15, 240, 51, 204, 60}, {0, 511, 511, 15, 240, 305, 450, 326}});

double qsqrt(double E) {
    E = sqrt(E);
    if (!params.quantize) {
        return E;
    }
    return floor(E*params.quantize_count+0.5)/params.quantize_count;
}

/* Numerical derivative (class double) */
/*
double fit_error_deriv(const Array<double> &in, Array<double> &out, shared_ptr<Filter<double> > f, const Array<double> &target, vector<double> &deriv) {
    double f0 = fit_error<double>(in, out, f, target);
    vector<double *> coeffs = f->coeffs();
    
    double eps = 1e-8;
    
    deriv.resize(coeffs.size());
    for (int i = 0; i < (int) coeffs.size(); i++) {
        double h = eps * (abs(*coeffs[i])+1);
        double orig = *coeffs[i];
        *coeffs[i] += h;
        
        double fp = fit_error<double>(in, out, f, target);
        *coeffs[i] = orig;
        deriv[i] = (fp - f0) / h;
    }
    return f0;
}
*/

void Profiling::init() {
    iter_changed = 0;
    current_iter = 0;
    max_iterations = 1;
    total_of_calls = 0;
    adaptive_iter.clear();
    mask_samples.clear();
    skip_count = 0;
}

void Profiling::reset() {
    T_fit = 0.0;
    T_bfgs = 0.0;
    T_deriv = 0.0;
    of_calls = 0;
    bfgs_calls = 0;
    bfgs_iters = 0;
    num_solve_initial = 0;
    num_solve_add = 0;
    add_count = 0;
}
    
Profiling profiling;

/* Automatic differentiation (class adouble) */
double fit_error_deriv(const Array<adouble> &in, Array<adouble> &out, shared_ptr<Filter<adouble> > f, const Array<double> &target, vector<double> &deriv) {
    profiling.of_calls++;
    profiling.total_of_calls++;
    vector<adouble *> coeffs = f->coeffs();
    deriv.resize(coeffs.size());
    
    stack.new_recording();
    //if (!adept::active_stack()->gradients_are_initialized()) {
    //    printf("gradients_are_initialized: %d\n", adept::active_stack()->gradients_are_initialized());
    //}
    //std::cout << stack;
    double T0_fit = wall_time();
    adouble y = fit_error<adouble>(in, out, f, target);
    profiling.T_fit += wall_time()-T0_fit;
    double T0_deriv = wall_time();
    y.set_gradient(1.0);
    //std::cout << stack;
    stack.compute_adjoint();
    
    for (int i = 0; i < (int) coeffs.size(); i++) {
        deriv[i] = coeffs[i]->get_gradient();
    }
    double ans = y.value();
    profiling.T_deriv += wall_time()-T0_deriv;
    return ans;
}

int translate_index(int i, const vector<MaskType> &mask, int dx, int dy) {
    int n = mask.size()/2;
    int row = int(floor(sqrt(n)));
    int imod = i % n;
    int shift = dx + dy*row;
    imod += shift;
    if (imod >= 0 && imod < n) {
        return i + shift;
    }
    return -1;
}

#define ABS(x) ((x)<0?(-(x)):(x))

vector<int> solve_translate(const vector<int> &dx_minL, const vector<int> &dx_maxL, const boxes_t &boxes) {
//    ASSERT(typeid(*f->filter) == typeid(FilterDAG<real>), "translate_border expected FilterSparse(FilterDAG())");
//    shared_ptr<FilterDAG<real> > dag = dynamic_pointer_cast<FilterDAG<real> >(f->filter);

//    vector<real *> coeffs(f->filter->coeffs());
//    ASSERT(dag->L.size() == 2, "expected 2 size DAG");
    ASSERT(boxes.size() == 2, "expected 2 boxes");
    ASSERT(dx_minL.size() == boxes.size(), "expected dx_minL size == boxes size");
    ASSERT(dx_maxL.size() == boxes.size(), "expected dx_maxL size == boxes size");
    
    vector<int> ans;
    int max_dx_abs = MAX(ABS(dx_minL[0]), ABS(dx_maxL[0]));
    for (int dx_abs = 0; dx_abs <= max_dx_abs; dx_abs++) {
        for (int dir = -1; dir <= 1; dir+=2) {
            int dx0 = dx_abs * dir;
            if (dx0 < dx_minL[0] || dx0 > dx_maxL[0]) { continue; }
//    for (int dx0 = dx_minL[0]; dx0 <= dx_maxL[0]; dx0++) {
            int dx1 = -dx0;
            if (dx1 >= dx_minL[1] && dx1 <= dx_maxL[1]) {
                if (params.verbose >= 2) {
                    printf("solve_translate returning %d %d\n", dx0, dx1);
                }
                ans.push_back(dx0);
                ans.push_back(dx1);
                return ans;
            }
        }
    }
    
    if (params.verbose >= 2) {
        printf("solve_translate no solution\n");
    }
    ans.push_back(0);
    ans.push_back(0);
    return ans;
}

int count_zero_row_col(shared_ptr<Array<int> > box, const vector<MaskType> &mask, int x0, int y0, int span_dx, int span_dy, int slab_dx, int slab_dy) {
    int ans = 0;
    int w = box->width();
    int h = box->height();
    if ((unsigned) x0 >= (unsigned) w || (unsigned) y0 >= (unsigned) h) { return 0; }
    for (;;) {
        int x = x0, y = y0;
        bool success = true;
        for (;;) {
            if ((unsigned) x >= (unsigned) w || (unsigned) y >= (unsigned) h) {
                fprintf(stderr, "%d, %d (x0,y0=%d,%d) out of bounds %d, %d\n", x, y, x0, y0, w, h);
                ASSERT(false, "count_zero_row_col out of bounds");
            }
            int idx = (*box)(y, x);
            ASSERT((unsigned) idx < (unsigned) mask.size(), "idx out of bounds");
            MaskType m = mask[idx];
            if (m == MASK_VAR) { success = false; break; }
            x += span_dx;
            y += span_dy;
            if ((unsigned) x >= (unsigned) w || (unsigned) y >= (unsigned) h) { break; }
        }
        if (!success) { break; }
        ans++;
        x0 += slab_dx;
        y0 += slab_dy;
        if ((unsigned) x0 >= (unsigned) w || (unsigned) y0 >= (unsigned) h) { break; }
    }
    return ans;
}

bool has_symmetry(const Array<double> &A, int dx, int dy, double epsilon=1e-5) {
    for (int y = 0; y < A.height(); y++) {
        for (int x = 0; x < A.width(); x++) {
            int xp = x, yp = y;
            if (dx) { xp = A.width()-1-x; }
            if (dy) { yp = A.height()-1-y; }
            double a = A(y, x);
            double b = A(yp, xp);
            if (fabs(a-b) > epsilon) {
                return false;
            }
        }
    }
    return true;
}

bool is_symmetry() {
    return params.allow_symmetry && (params.symmetry_h || params.symmetry_v);
}

void detect_symmetry(const Array<double> &target) {
    if (params.symmetry) {
        params.symmetry_h = has_symmetry(target, 1, 0);
        params.symmetry_v = has_symmetry(target, 0, 1);
        printf("Symmetry detection: horizontal: %d, vertical: %d\n", int(params.symmetry_h), int(params.symmetry_v));
    } else {
        params.symmetry_v = params.symmetry_h = false;
    }
//    exit(1);
}

int sample_box(boxes_t boxes) {
    vector<int> box_idxL;
    for (int i = 0; i < (int) boxes.size(); i++) {
        for (int j = 0; j < boxes[i]->nelems; j++) {
            box_idxL.push_back(i);
        }
    }
    if (box_idxL.size() == 0) { return -1; }
    //ASSERT(box_idxL.size(), "boxes should be non-empty");
    return box_idxL[rand()%box_idxL.size()];
}

vector<int> get_multisize_sizes() {
    if (!params.multisize || params.multisize_sizes <= 1) {
        return {params.filter_w};
    } else {
        vector<int> ans;
        double min_size = 3;
        double max_size = params.filter_w;
        double scale = max_size / min_size;
        for (int i = 0; i < params.multisize_sizes; i++) {
            double t = i * 1.0 / (params.multisize_sizes-1);
            double sz = min_size * pow(scale, t) + 1e-6;
            int szi = int(sz);
            while (szi % 2 == 0 || (ans.size() > 0 && szi <= ans[ans.size()-1])) { szi++; }
            ans.push_back(szi);
        }
        return ans;
    }
}

bool terminate_early(const vector<double> &Tavg_L, double T_start_optim) {
    if (params.terminate_early && (int) Tavg_L.size() >= params.terminate_iters) {
        double Tavg_max = 0.0;
        for (int i = 0; i < params.terminate_iters; i++) {
            Tavg_max = MAX(Tavg_max, Tavg_L[Tavg_L.size()-1-i]);
        }
        double ratio = Tavg_L[Tavg_L.size()-1] / Tavg_max;
        if (ratio >= 1-params.terminate_epsilon) {
            return true;
        }
    }
    if (params.max_time > 0 && (wall_time()-T_start_optim)/60.0 > params.max_time) {
        return true;
    }
    return false;
}

string mask_to_str(const vector<MaskType> &m) {
    string ans(m.size(), '0');
    for (int i = 0; i < (int) m.size(); i++) {
        if (m[i] == MASK_VAR) { ans[i] = '1'; }
    }
    return ans;
}

/*
Array<double> spatial_weights_w;
bool spatial_weights_init = false;

void get_spatial_weights(const Array<double> &target, Array<double> *&w_ptr) {
    if (!spatial_weights_init) {
        spatial_weights_w.resize(target.sizes);
        int height = spatial_weights_w.height(), width = spatial_weights_w.width();
        double sigma = MAX(width, height) * params.spatial_sigma_frac;
        double sigma2 = sigma * sigma;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int dx = x - width/2;
                int dy = y - height/2;
                spatial_weights_w(y,x) = 1-exp(-(dx*dx+dy*dy)/(2.0*sigma2));
            }
        }
        spatial_weights_w = normalized(spatial_weights_w);
        spatial_weights_init = true;
    }
    ASSERT(spatial_weights_w.sizes == target.sizes, "w size != target size");
    w_ptr = &spatial_weights_w;
}
*/