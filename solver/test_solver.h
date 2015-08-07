
#ifndef _test_solver
#define _test_solver

#include "solver.h"
#include "simplify.h"
#include "mutate.h"
#include "filter_json.h"
#include <sys/stat.h>

void test_filters();

#define DEFAULT_TARGET_W  19        /* Was 15 */
#define REPEATS           10000

template<class real>
double test_fit() {
    int w = 7;
    int Kw = 5;
    srand(0);
    
    double sigma = 1.0;
    Array<real> in = impulse<real>(w, w);
    auto out = in;
    Array<double> target({w, w});
    for (int y = 0; y < w; y++) {
        for (int x = 0; x < w; x++) {
            int dx = x-w/2;
            int dy = y-w/2;
            target(y,x) = exp(-(dx*dx+dy*dy)/(2*sigma*sigma));
        }
    }
    
    auto KA = Array<double>::random({Kw, Kw});
    auto KB = Array<double>::random({Kw, Kw});
    
    FilterFIR<real> fA(KA);
    FilterFIR<real> fB(KB);
    FilterDAG<real> f({&fA, &fB});

    vector<double> deriv;
    
    double Tval = 0;
    
    for (int k = 0; k < 2; k++) {
        double T0 = wall_time();
        double fval = 0;
        for (int iter = 0; iter < REPEATS; iter++) {
            fval = fit_error_deriv(in, out, &f, target, deriv);
//            fval = fit_error_deriv(in, out, &fA, target, deriv);
            //printf("%f\n", to_double(err));
        }
        Tval = (wall_time() - T0)/REPEATS;

        printf("Time: %e secs\n", Tval);
        printf("f = %f\n", fval);

        for (int i = 0; i < (int) deriv.size(); i++) {
            printf("%f ", deriv[i]);
        }
        printf("\n\n");
    }
    return Tval;
}

#define TEST_PREAMBLE() \
int w = DEFAULT_TARGET_W; \
int Kw = 7; \
Array<real> in(vector<int>{w, w}); \
in.clear(); \
in(w/2,w/2) = 1; \
Array<real> out(vector<int>{w, w}); \
Array<double> target(vector<int>{w, w}); \
shared_ptr<Filter<real> > f;

#define GET_TEST_KERNELS() \
auto a = shared_ptr<FilterFIR<real> >(new FilterFIR<real>(Array<double>::random({Kw, Kw}))); \
auto b = shared_ptr<FilterFIR<real> >(new FilterFIR<real>(Array<double>::random({Kw, Kw}))); \
auto a_filt = static_pointer_cast<Filter<real> >(a); \
auto b_filt = static_pointer_cast<Filter<real> >(b); \
\
double sigma = 2; \
for (int y = 0; y < w; y++) { \
    for (int x = 0; x < w; x++) { \
        int dx = x-w/2; \
        int dy = y-w/2; \
        target(y,x) = exp((-dx*dx-dy*dy)/(2*sigma*sigma)); \
    } \
} \
target = normalize_target(target); \
\
f = shared_ptr<Filter<real> >(new FilterDAG<real>({a_filt, b_filt}));

#define GET_TEST_KERNELS2() \
auto a = shared_ptr<FilterFIR<real> >(new FilterFIR<real>(Array<double>::random({Kw, Kw}))); \
auto b = shared_ptr<FilterIIR<real> >(new FilterIIR<real>(zeros<double>({Kw}), zeros<double>({Kw}))); \
b->K(0) = 1.0; \
b->F(0) = 0.5; \
auto a_filt = static_pointer_cast<Filter<real> >(a); \
auto b_filt = static_pointer_cast<Filter<real> >(b); \
\
for (int y = 0; y < w; y++) { \
for (int x = 0; x < w; x++) { \
int dx = x-w/2; \
int dy = y-w/2; \
target(y,x) = abs(dx) <= 3 && abs(dy) <= 3; \
} \
} \
target = normalize_target(target); \
\
f = shared_ptr<Filter<real> >(new FilterDAG<real>({a_filt, b_filt}));


#define RANDOMIZE_TEST_KERNELS() \
a->K = Array<double>::random({Kw, Kw}); \
b->K = Array<double>::random({Kw, Kw});


template<class real>
void test_fit_bfgs() {
    vector<Array<real> > temp_images;

    srand(0); //srand(time(0));

    TEST_PREAMBLE();
    GET_TEST_KERNELS();
    //f->apply(in, out);
    //printf("f applied in:\n%s\n\n", in.str().c_str());
    //printf("f applied out:\n%s\n\n", out.str().c_str());

    vector<double> deriv;
    //out.clear();
    double Ta = wall_time();
    double finitial = to_double(fit_error<real>(in, out, f, target));
    double Tb = wall_time();
    //printf("f applied out (via fit_error_deriv):\n%s\n\n", out.str().c_str());
    printf("f initial: %f\n", finitial);

    //printf("d\n"); fflush(stdout);
    int iters = 0;
    double T0 = wall_time();
    double ffinal = optimize_filter_bfgs(in, out, f, target, &iters);
    double T1 = wall_time();
    f->apply({&in}, out, temp_images);
    double T2 = wall_time();
    printf("\noutput: %s\n", normalize_target(out).str().c_str());
    printf("\ntarget: %s\n", normalize_target(target).str().c_str());
    printf("\nOptimized in %d iters, %e secs, err=%f\n", iters, T1-T0, ffinal);
    printf("apply time: %e secs\n", T2-T1);
    printf("fit_error time (adouble or double): %e secs\n", Tb-Ta);
    
    /*
    auto a_double = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(a->K));
    auto b_double = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(b->K));
    auto f_double = shared_ptr<FilterDAG<double> >(new FilterDAG<double>({a_double, b_double}));
    auto f_base = static_pointer_cast<Filter<double> >(f_double);
    auto in_double = Array<double>(in);
    auto out_double = Array<double>(in);
    
    double T3 = wall_time();
    //fit_error<double>(in_double, out_double, f_base, target);
    double T4 = wall_time();
    printf("fit_error time (double): %e secs\n", T4-T3);
    printf("\n");
    */
}

template<class real>
double get_1norm_scale(const Array<real> &in, Array<real> &out, shared_ptr<FilterSparse<real> > filter, const Array<double> &target) {
    fit_error<real>(in, out, filter, target, false);     // Do not rescale output
    double sum = 0.0;
    for (int y = 0; y < (int) out.height(); y++) {
        for (int x = 0; x < (int) out.width(); x++) {
            sum += fabs(out(y, x).value());
        }
    }
    return sum > 0 ? 1.0/sum: 0.0;
}

template<class real>
void write_pareto(const Array<real> &in, Array<real> &out, const Array<double> &target, PARETO_TYPE(real) *pareto, bool print=true, bool check=true, string solver_cmd="", bool scale_1norm=false) {
    if (print) {
        printf("Pareto:\n%s\n", pareto->str().c_str());
    }
    
    if (params.pareto.size()) {
        FILE *f = fopen(params.pareto.c_str(), "wt");
        fprintf(f, pareto->str().c_str());
        fclose(f);
    }
    
    if (params.pareto_full.size()) {
        FILE *f = fopen(params.pareto_full.c_str(), "wt");
        fprintf(f, "[\n");
        for (int i = 0; i < (int) pareto->L.size(); i++) {
            auto filter = pareto->L[i].p;
            filter->zero_masked_coeffs();
            
            if (params.flatten) {
                ASSERT2(typeid(*filter) == typeid(FilterSparse<real>), "expected filter to be FilterSparse");
                shared_ptr<FilterDAG<real> > dag = static_pointer_cast<FilterDAG<real> >(filter->filter);
                dag = dag->flatten();
                filter = make_shared<FilterSparse<real> >(dag);
                filter->mask_from_coeffs();
            }
            
//            FILE *f0 = fopen("f0.txt", "at");
//            fprintf(f0, "Writing Pareto to %s, index %d\n", params.pareto_full.c_str(), i);
//            fclose(f0);
            
            string s = filter->str(false);
            s = str_replace(s, "\n", "\n    ");
            
            string mask_str = vector_to_str_int<MaskType>(filter->mask);
            real filter_scale = 0.0;
            vector<double> features = filter->features();
            
            RECLAIM_ADEPT_MEMORY();
            double err_precise = sqrt(to_double(fit_error<real>(in, out, filter, target, true, &filter_scale)));     // Rescale output
            if (scale_1norm || params.scale_1norm) {
                filter_scale = get_1norm_scale(in, out, filter, target);
            } else {
                filter_scale *= params.target_orig_norm;
            }
            fprintf(f, "{\n"
                    "    \"name\": \"program%03d\",\n"
                    "    \"time\": %.17g,\n"
                    "    \"error\": %.17g,\n"
                    "    \"mask\": %s,\n"
                    "    \"solver_cmd\": \"%s\",\n"
                    "    \"features\": %s,\n"
                    "    \"program\": %s,\n"
                    "    \"scale\": %.16g,\n"
                    "    \"out\": %s\n"
                    "}", i, pareto->L[i].T, err_precise, mask_str.c_str(), solver_cmd.c_str(), vector_to_str_real(features).c_str(), s.c_str(), to_double(filter_scale), out.str().c_str());
            if (i < (int) (pareto->L.size()-1)) {
                fprintf(f, ",");
            }
            fprintf(f, "\n");
        }
        fprintf(f, "]\n");
        fclose(f);
    }
    if (check) {
        check_pareto(in, out, target, pareto);
    }
}

template<class real>
void test_fit_sparse() {
    TEST_PREAMBLE();
    GET_TEST_KERNELS();

    if (params.kernel_function == 1) {
        double sigma = 2;
        for (int y = 0; y < w; y++) {
            for (int x = 0; x < w; x++) {
                int dx = x-w/2;
                int dy = y-w/2;
                target(y,x) = exp((-dx*dx-dy*dy)/(2*sigma*sigma))*cos(dx/4.0*dy/4.0);
            }
        }
    } else if (params.kernel_function == 2) {
/*        for (int y = 0; y < w; y++) {
            for (int x = 0; x < w; x++) {
                int dx = x-w/2;
                int dy = y-w/2;
                target(y,x) = (abs(dx) <= 5 && abs(dy) <= 5) ? 1: 0;
            }
        }*/
        GET_TEST_KERNELS2();
    }
    
    auto pareto = new PARETO_TYPE(real)();
    auto fsparse = shared_ptr<FilterSparse<real> >(new FilterSparse<real>(f));
    
    if (params.init_correct_sparsity) {
        for (int i = 0; i < fsparse->mask.size(); i++) {
            fsparse->mask[i] = MASK_ZERO;
        }
        for (int dy = 0; dy < Kw; dy++) {
            fsparse->mask[dy*Kw+Kw/2] = MASK_VAR;
            fsparse->mask[Kw*Kw+dy+(Kw/2)*Kw] = MASK_VAR;
        }
    }
    
    for (int iter = 0; iter < params.random_iters; iter++) {
        RANDOMIZE_TEST_KERNELS();
        optimize_sparse_filter(in, out, fsparse, target, pareto);
    }
    //optimize_sparse_filter_expand(in, out, fsparse, target, pareto);
    
    write_pareto(in, out, target, pareto);
}

void test_solver();

#define FILTER_MAP(real) map<int, vector<shared_ptr<FilterDAG<real> > > >

template<class real>
void get_combinatoric_filters(FILTER_MAP(real) &filters, bool check_subprograms=false) {
    auto seen = new set<string>();
    (filters)[0] = {shared_ptr<FilterDAG<real> >(new FilterDAG<real>(vector<shared_ptr<Filter<real> > >({})))};
    for (int count = 0; count < params.max_combinatoric; count++) {
        for (int i = 0; i < (int) (filters)[count].size(); i++) {
            vector<shared_ptr<FilterDAG<real> > > addL = add_mutate((filters)[count][i], *seen, true);
            //printf("addL size: %d\n", int(addL.size()));
            for (int j = 0; j < (int) addL.size(); j++) {
                int c = addL[j]->L.size();
                (filters)[c].push_back(addL[j]);
                if (check_subprograms) {
                    get_subprogram_fragment(addL[j]);
                }
            }
        }
    }
    delete seen;
}

template<class real>
void get_fixed_topology(FILTER_MAP(real) &filters) {
    ASSERT(filters.size() == 0, "-topology expects filters map to initially be size 0");
    
    shared_ptr<FilterDAG<real> > f;
    if (params.fixed_topology == "fir_fir") {
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto b = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "du") {
        auto a = make_shared<FilterDownsample<real> >();
        auto b = make_shared<FilterUpsample<real> >();
        vector<shared_ptr<Filter<real> > > fL({a, b});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "d_fir_u") {
        auto a = make_shared<FilterDownsample<real> >();
        auto b = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto c = make_shared<FilterUpsample<real> >();
        vector<shared_ptr<Filter<real> > > fL({a, b, c});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "fir_d_fir_u_fir") {
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto b = make_shared<FilterDownsample<real> >();
        auto c = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto d = make_shared<FilterUpsample<real> >();
        auto e = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b, c, d, e});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "twoband") {
        set<string> seen;
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        vector<shared_ptr<Filter<real> > > fL({a});
        auto f0 = make_shared<FilterDAG<real> >(fL);
        f = add_mutate_band(f0, seen, true)[0];
    } else if (params.fixed_topology == "threeband") {
        set<string> seen;
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        vector<shared_ptr<Filter<real> > > fL({a});
        auto f0 = make_shared<FilterDAG<real> >(fL);
        auto f1 = add_mutate_band(f0, seen, true)[0];
        f = add_mutate_band(f1, seen, true)[0];
    } else if (params.fixed_topology == "fir_d_fir_u") {
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto b = make_shared<FilterDownsample<real> >();
        auto c = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto d = make_shared<FilterUpsample<real> >();
        vector<shared_ptr<Filter<real> > > fL({a, b, c, d});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "firh_firv") {
        auto a = make_shared<FilterFIR<real> >(fir_initial_params_h<real>());
        auto b = make_shared<FilterFIR<real> >(fir_initial_params_v<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "fir_iir") {
        auto a = make_shared<FilterFIR<real> >(fir_initial_params<real>());
        auto b = make_shared<FilterIIR<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "iir_iir") {
        auto a = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "iir4") {
        auto a = make_shared<FilterIIR<real> >      (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b = make_shared<FilterIIRMinusX<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto c = make_shared<FilterIIRPlusY<real> > (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto d = make_shared<FilterIIRMinusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        vector<shared_ptr<Filter<real> > > fL({a, b, c, d});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "iir_iir2") {
        auto a1 = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b1 = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto a2 = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b2 = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        vector<shared_ptr<Filter<real> > > fL({a1, b1, a2, b2});
        f = make_shared<FilterDAG<real> >(fL);
    } else if (params.fixed_topology == "iir_iir3") {
        auto a1 = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b1 = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto a2 = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b2 = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto a3 = make_shared<FilterIIR<real> >     (iir_initial_params_K<real>(), iir_initial_params_F<real>());
        auto b3 = make_shared<FilterIIRPlusY<real> >(iir_initial_params_K<real>(), iir_initial_params_F<real>());
        vector<shared_ptr<Filter<real> > > fL({a1, b1, a2, b2, a3, b3});
        f = make_shared<FilterDAG<real> >(fL);
    } else {
        fprintf(stderr, "Unknown topology %s\n", params.fixed_topology.c_str());
        exit(1);
    }
    if (params.verbose) {
        printf("Initial DAG:\n%s\n\n", f->str().c_str());
    }
    (filters)[f->L.size()].push_back(f);
}

template<class real>
int filter_count(FILTER_MAP(real) &filters, bool verbose=true) {
    if (verbose) {
        printf("Filter count:\n");
    }
    int ans = 0;
    vector<double> sizes;
    for (auto it = filters.begin(); it != filters.end(); it++) {
        if (verbose) {
            int n_twolevel = 0;
            int n_threelevel = 0;
            for (auto filter: it->second) {
                if (filter->L.size()) {
                    int levels = dag_levels(filter, sizes);
                    if (levels == 1) { n_twolevel++; }
                    else if (levels == 2) { n_threelevel++; }
                }
            }
            printf("  Filters size %d: %d (two level: %d, three level: %d)\n", it->first, (int) it->second.size(), n_twolevel, n_threelevel); //(int) (*filters)[count].size());
            ans += (int) it->second.size();
        }
    }
    printf("\n");
    return ans;
}

vector<string> parse_args(int argc, char *argv[]);
void set_out_prefix(const string &a);

template<class real>
void test_random_mutate(int nfilters=200, int nsteps=5) {
    vector<shared_ptr<FilterDAG<real> > > prevL;
    int crossover_count = 0;
    int mutate_band_count = 0;
    int delete_mutate_count = 0;
    int add_mutate_count = 0;
    int replace_mutate_count = 0;
    for (int i = 0; i < nfilters; i++) {
        if (params.verbose) { printf("%d/%d\n", i, nfilters); }
        set<string> seen;
        shared_ptr<FilterDAG<real> > f;
        if (rand() % 2 == 0) {
            f = make_shared<FilterDAG<real>>(vector<shared_ptr<Filter<real> > >({}));
        } else {
            vector<string> topology_strL = {"fir_fir", "firh_firv", "fir_iir", "iir_iir", "iir4", "iir_iir2", "iir_iir3"};
            params.fixed_topology = topology_strL[rand()%int(topology_strL.size())];
            FILTER_MAP(real) fmap;
            get_fixed_topology(fmap);
            bool success = false;
            for (auto p: fmap) {
                if (p.second.size()) {
                    f = p.second[0];
                    success = true;
                }
            }
            ASSERT2(success, "Failed to get fixed topology");
        }
        for (int j = 0; j < nsteps; j++) {
            get_subprogram_fragment(f);
            vector<shared_ptr<FilterDAG<real> > > fL;
            seen.clear();
            while (1) {
                int r = rand() % 3;
                if (r == 0 && f->L.size() > 2) {
                    fL = delete_mutate<real>(f, seen, false);
                    delete_mutate_count++;
                } else if (r == 1 && f->L.size() >= 1) {
                    fL = replace_mutate<real>(f, seen);
                    replace_mutate_count++;
                } else {
                    fL = add_mutate<real>(f, seen, false);
                    add_mutate_count++;
                }
                if (fL.size()) { break; }
            }
            //if (!fL.size()) {
            //    printf("Could not mutate original, using r=%d:\n%s\n\n", r, f->str().c_str());
            //}
            //ASSERT2(fL.size(), "Failed to mutate");
            ASSERT2(fL.size() == 1, "Got more than 1 mutation");
            f = fL[0];
            if (prevL.size() == 0 || f->L.size() <= 20) {
                prevL.push_back(f);
            }
            fL = crossover(f, prevL[rand()%int(prevL.size())], seen);
            if (fL.size()) {
                crossover_count++;
                f = fL[0];
            }
            fL = add_mutate_band(f, seen);
            if (fL.size()) {
                mutate_band_count++;
                f = fL[0];
            }
        }
    }
    ASSERT2(delete_mutate_count, "no delete_mutate");
    ASSERT2(add_mutate_count, "no add mutate");
    ASSERT2(replace_mutate_count, "no replace mutate");
    ASSERT2(crossover_count, "no crossover");
    ASSERT2(mutate_band_count, "no add_mutate_band");
    printf("add_mutate/delete_mutate(all=false):   OK\n  (%d crossover, %d add_mutate_band, %d delete_mutate, %d add_mutate, %d replace_mutate)\n", crossover_count, mutate_band_count, delete_mutate_count, add_mutate_count, replace_mutate_count);
}

template<class real>
void test_mutate(int argc, char *argv[]) {
    parse_args(argc, argv);

    test_random_mutate<real>();

    double T0 = wall_time();
    FILTER_MAP(real) filters;
    Array<double> target(vector<int>({10, 10}));
    
    get_combinatoric_filters(filters, true);
    //for (int count = 0; count <= MAX_COUNT; count++) {
    filter_count(filters);
    //printf("dealloc seen\n");
    //printf("done\n");
    //printf("dealloc filters\n");
    printf("add_mutate(all=true) in %f secs\n", wall_time()-T0);
    //printf("done\n");
}

void usage();

//Pareto<shared_ptr<FilterSparse<real> > >
template<class real>
void set_pareto_level(PARETO_TYPE(real) *pareto, int level, const Array<real> &in, Array<real> &out, const Array<double> &target) {
    check_pareto(in, out, target, pareto);
    double eps = 1e-6;
    for (int i = 0; i < (int) pareto->L.size(); i++) {
        auto fval = dynamic_pointer_cast<Filter<real> >(pareto->L[i].p);
        double err0 = fit_error_double_sqrt(in, out, fval, target);

        pareto->L[i].p->set_level(level);
        double T1 = pareto->L[i].T;
        double T2 = pareto->L[i].p->time();
        if (fabs(T1 - T2) >= eps) {
            fprintf(stderr, "times differ in set_pareto_level, %d/%d: %f %f\n", i, (int) pareto->L.size(), T1, T2);
            ASSERT(0, "times differ in set_pareto_level");
        }
        double err1 = pareto->L[i].err;
        double err2 = fit_error_double_sqrt(in, out, fval, target);
        if (fabs(err1 - err2) >= eps) {
            fprintf(stderr, "err differ in set_pareto_level, %d/%d: err0=%f err1=%f err2=%f, T1=%f T2=%f\n", i, (int) pareto->L.size(), err0, err1, err2, T1, T2);
            ASSERT(0, "err differ in set_pareto_level");
        }
    }
}

template<class real>
void convolve(const Array<real> &Aorig0, const Array<real> &K, Array<real> &Aout) {
    Array<real> Aorig(Aorig0);
    FilterFIR<real> f(K);
    vector<Array<real> > temp_images;
    
    f.apply({&Aorig}, Aout, temp_images);
}

template<class real>
void get_target_in_out(string target_filename, Array<double> &target, Array<real> &in, Array<real> &out) {
    target = load_matrix<double>(target_filename);
    params.target_orig_norm = norm_target(target);
    ASSERT(target.sizes.size() == 2, "expected 2D target");
    target = normalize_target(target);
    if (params.in.size()) {
        in = load_matrix<double>(params.in);
    } else {
        in = impulse<real>(target.width(), target.height());
    }
    out = in;
    
    if (params.target_w > 0) {
        target = zero_pad_center(target, params.target_w, params.target_w);
        in     = zero_pad_center(in,     params.target_w, params.target_w);
        out    = zero_pad_center(out,    params.target_w, params.target_w);
    }
    
    if (params.preconvolve.size()) {
        printf("Preconvolving with %s\n", params.preconvolve.c_str());
        Array<double> P  = load_matrix<double>(params.preconvolve);
        Array<real> Pr = load_matrix<real>(params.preconvolve);
        
        convolve(in, Pr, in);
        convolve(target, P, target);
    }
}

#define GET_TARGET_IN_OUT(target_filename) \
Array<double> target; \
Array<real> in, out; \
get_target_in_out<real>(target_filename, target, in, out);

template<class real>
shared_ptr<PARETO_TYPE(real)> load_pareto(string json_filename, const Array<real> &in, Array<real> &out, const Array<double> &target, bool recalc_solve=false, bool retain_all=false) {
    vector<shared_ptr<FilterSparse<real> > > L = json_to_filters<real>(json_filename);
    auto pareto_main = make_shared<PARETO_TYPE(real)>();
    for (int i = 0; i < (int) L.size(); i++) {
        auto f = dynamic_pointer_cast<Filter<real> >(L[i]);
        auto fsparse(L[i]);
//        printf("original mask: %s\n", vector_to_str_int(fsparse->mask).c_str());
        double err = fit_error_double_sqrt(in, out, f, target);
        double time = f->time();
//        auto fsparse = make_shared<FilterSparse<real> >(f);
//        fsparse->mask_from_coeffs();
        if (recalc_solve) {
            vector<double> gradient;
            optimize_and_add(in, out, fsparse->copy_sameclass(), target, gradient, &(*pareto_main), NULL, 0, SPARSE_NULL);
        } else {
            if (!params.retain_all && !retain_all) {
                pareto_main->add(time, err, fsparse);
            } else {
                pareto_main->L.push_back(ParetoPoint<shared_ptr<FilterSparse<real> > >(time, err, fsparse));
            }
        }
    }
    return pareto_main;
}

template<class real>
void main_recalc(int argc, char *argv[]) {
    if (argc < 2) {
        usage();
    }
    string target_filename = string(argv[0]);
    string json_filename = string(argv[1]);
    parse_args(argc-2, argv+2);
    
    GET_TARGET_IN_OUT(target_filename);
    
    auto pareto_main = load_pareto(json_filename, in, out, target, params.recalc_solve);
    write_pareto(in, out, target, &(*pareto_main));
}

void print_begin_solver(int argc, char *argv[]);

template<class real>
class FilterContext { public:
    const Array<real> in;
    Array<real> out;
    const Array<real> target;
    
    int nfilters;
    double best_max_time;
    int argc;
    char **argv;
    string target_filename;

    FilterContext(const Array<real> &in_, const Array<real> &out_, const Array<real> &target_, int argc_, char **argv_, string target_filename_)
        :in(in_), out(out_), target(target_), nfilters(0), best_max_time(1e100), argc(argc_), argv(argv_), target_filename(target_filename_) { }
};

#define FILTER_CONTEXT shared_ptr<FilterContext<real> >

template<class real>
double filter_min_time(shared_ptr<FilterDAG<real> > f0) {
    auto fcopy = f0->copy_sameclass();
    vector<real *> f_coeffs = fcopy->coeffs();
    for (int i = 0; i < (int) f_coeffs.size(); i++) {
        *f_coeffs[i] = 0.0;
    }
    return fcopy->time();
}

template<class real>
class GAIndividual { public:
    shared_ptr<PARETO_TYPE(real)> pareto;
    shared_ptr<FilterDAG<real> > f;             /* Initial FilterDAG with the right topology (before solver) */
    string simplify_str;
    double avg_time;
    double max_time;
    double min_time;
    double objective;           /* Fitness objective (to be minimized) */
    string out_path;
    FILTER_CONTEXT context;
    
    void after_set_pareto() {
        //write_pareto(context->in, context->out, context->target, &(*pareto), false);
        ASSERT2(pareto->L.size(), "pareto size should be nonzero");
        avg_time = pareto->avg_time(params.pareto_error_thresh);
        max_time = pareto->max_time(params.pareto_error_thresh);
        auto fsparse = pareto->L[0].p;
        ASSERT2(typeid(*fsparse->filter) == typeid(FilterDAG<real>), "expected Pareto to contain FilterSparse(FilterDAG())");
        auto f = dynamic_pointer_cast<FilterDAG<real> >(fsparse->filter);
        min_time = filter_min_time(f);
        context->best_max_time = MIN(context->best_max_time, max_time);
    }
    
    GAIndividual(FILTER_CONTEXT context_, shared_ptr<FilterDAG<real> > f_)
        :f(f_), context(context_) {
        simplify_str = simplify(f);
        objective = 0.0;
//        get_pareto();
    }
};

template<class real>
class GAIndividualCompare { public:
    bool operator ()(const shared_ptr<GAIndividual<real> > &a, const shared_ptr<GAIndividual<real> > &b) {
        return a->objective < b->objective;
    }
};

#define CHECK_IIR_CLASS(cls) \
        if (typeid(*f0->L[i].f) == typeid(FilterIIR<real>)) { \
            auto fcurrent = static_pointer_cast<FilterIIR<real> >(f0->L[i].f); \
            fcurrent->K = iir_initial_params_K<real>(wp); \
            fcurrent->F = iir_initial_params_F<real>(wp); \
        }

#define VERBOSE_MSG(s) if (params.verbose) { printf(s); }
template<class real>
bool initial_check(FILTER_CONTEXT context, shared_ptr<FilterDAG<real> > f0, bool check_time=true) {
    if (params.verbose) { printf("initial_check:\n%s\n", f0->str(false).c_str()); }
    if (f0->L.size() == 0) { VERBOSE_MSG("initial_check failed: size is zero\n"); return false; }
    if (f0->coeffs().size() == 0) { VERBOSE_MSG("initial_check failed: coeffs is empty\n"); return false; }
    
    vector<double> sizes;
    int levels = dag_levels(f0, sizes);                /* Check that (1<<levels) divides output width */
    int w = context->out.width();
    if (w % (1<<levels) != 0) { VERBOSE_MSG("initial_check failed: 1<<levels does not divide w\n"); return false; }
    if (params.ga_always_downsample && levels == 0) { VERBOSE_MSG("initial_check failed: levels is zero and always_downsample\n"); return false; }
    
    if (check_time) {
        //ASSERT(pop_prev.size(), "expected pop_prev to be non-empty");
        double min_time = filter_min_time(f0);
        if (min_time > context->best_max_time) {
            if (params.verbose) { printf("initial_check failed: filter time too high (best_max_time: %f, min_time: %f)\n", context->best_max_time, min_time); }
            return false;                       /* If our min time exceeds any existing individual's max time, quit early. */
        }
    }
    
    ASSERT(sizes.size() == f0->L.size(), "sizes size should be f0->L.size()");
    for (int i = 0; i < (int) f0->L.size(); i++) {
        double frac = sqrt(sizes[i]);
        int wp = int(params.filter_w * frac);
        if (wp < 1) { wp = 1; }
        if (typeid(*f0->L[i].f) == typeid(FilterFIR<real>)) {
            auto fcurrent = static_pointer_cast<FilterFIR<real> >(f0->L[i].f);
            fcurrent->K = fir_initial_params<real>(wp);
        }
        CHECK_IIR_CLASS(FilterIIR<real>);
        CHECK_IIR_CLASS(FilterIIRMinusX<real>);
        CHECK_IIR_CLASS(FilterIIRPlusY<real>);
        CHECK_IIR_CLASS(FilterIIRMinusY<real>);
    }
    
    VERBOSE_MSG("initial_check: success\n");
    return true;
}

template<class real>
void sort_population(vector<shared_ptr<GAIndividual<real> > > &pop) {
    vector<int> sub_idx;
    for (int i = 0; i < (int) pop.size(); i++) { sub_idx.push_back(i); }
    vector<double> OF = get_tournament_objective(pop, sub_idx, int(pop.size()));
    for (int i = 0; i < (int) pop.size(); i++) { pop[i]->objective = OF[i]; }

    sort(pop.begin(), pop.end(), GAIndividualCompare<real>());
}

void run_commands(const vector<string> &L);

/* Evaluate Pareto frontiers for population then sort population */
template<class real>
void set_population_pareto(FILTER_CONTEXT context, vector<shared_ptr<GAIndividual<real> > > &pop) {
    if (params.verbose) {
        print_generation(context, pop, -1, "(before computing Pareto)");
        fflush(stdout);
    }

    vector<string> cmdL;
    for (int i = 0; i < (int) pop.size(); i++) {
        if (params.verbose) {
            printf("  get_pareto[%d/%d], %s (making next generation)\n", i, int(pop.size()), pop[i]->simplify_str.c_str());
        }
        shared_ptr<FilterSparse<real> > f0 = make_shared<FilterSparse<real> >(pop[i]->f);
        
        char buf[1024];
        sprintf(buf, "f%05d", context->nfilters);
        context->nfilters++;
        pop[i]->out_path = params.out_arg + buf;
//        pop[i]->set_out_prefix(pop[i]->out_path);

        string in_pareto = pop[i]->out_path + "_in.json";
        
        char cmd_prefix[1024];
        sprintf(cmd_prefix, "./solve_kernel solve %s %s", context->target_filename.c_str(), in_pareto.c_str());
        
        string cmd = cmd_prefix;
        for (int j = 0; j < context->argc; j++) {
            cmd += string(" ") + context->argv[j];
        }
        char cmd_suffix[1024];
        sprintf(cmd_suffix, " -out %s -verbose 0 &> %s_out.txt", pop[i]->out_path.c_str(), pop[i]->out_path.c_str());
        cmd += cmd_suffix;
        
        cmdL.push_back(cmd);
        
        PARETO_TYPE(real) pareto1;
        pareto1.add(0.0, 0.0, f0);
        string pareto0 = params.pareto;
        string pareto_full0 = params.pareto_full;
        params.pareto_full = in_pareto;
        params.pareto = "";
        write_pareto(context->in, context->out, context->target, &pareto1, false, false);
        params.pareto = pareto0;
        params.pareto_full = pareto_full0;
    }
    
    fflush(stdout);
    run_commands(cmdL);
    
    auto pop0 = pop;
    pop.clear();
    
    for (int i = 0; i < (int) pop0.size(); i++) {
        string json_filename = pop0[i]->out_path + "_pareto_full.txt";
        bool exists = true;
        FILE *f = fopen(json_filename.c_str(), "rt");
        if (!f) { exists = false; }
        else { fclose(f); }
        if (exists) {
            pop.push_back(pop0[i]);
            int pop_idx = pop.size()-1;
            pop[pop_idx]->pareto = load_pareto(json_filename, context->in, context->out, context->target);
            pop[pop_idx]->after_set_pareto();
        } else {
            printf("Warning: could not read %s\n", json_filename.c_str());
        }
    }
    sort_population(pop);
}

template<class real>
void init_population(FILTER_CONTEXT context, set<string> &seen, vector<shared_ptr<GAIndividual<real> > > &pop, bool set_pareto) {
    FILTER_MAP(real) filters;
    params.max_combinatoric = 3;
    double verbose0 = params.verbose;
    params.verbose = 0;
    RECLAIM_ADEPT_MEMORY();
    get_combinatoric_filters<real>(filters);
    params.verbose = verbose0;

    if (params.seed_convpyr) {
        string seed_file = params.out_arg + string("convpyr.json");
        char buf[4096];
        sprintf(buf, "python topologies.py initial %d %s %d", context->out.width(), seed_file.c_str(), params.filter_w);
        system(buf);
        
        //shared_ptr<PARETO_TYPE(real)> load_pareto(string json_filename, const Array<real> &in, Array<real> &out, const Array<double> &target, bool recalc_solve=false, bool retain_all=false)
        shared_ptr<PARETO_TYPE(real)> seedL = load_pareto(seed_file, context->in, context->out, context->target, false, true);
        for (int i = 0; i < (int) seedL->L.size(); i++) {
            shared_ptr<FilterSparse<real> > fsparse = seedL->L[i].p;
            ASSERT2(typeid(*fsparse->filter) == typeid(FilterDAG<real>), "expected FilterSparse(FilterDAG)");
            shared_ptr<FilterDAG<real> > fdag = static_pointer_cast<FilterDAG<real> >(fsparse->filter);
            
            string simplify_str = simplify(fdag);
            seen.insert(simplify_str);
            if (pop.size() < params.ga_population) {
                pop.push_back(make_shared<GAIndividual<real> >(context, fdag));
            }
        }
    }
    
    for (auto it = filters.begin(); it != filters.end(); it++) {
        for (int i = 0; i < (int) it->second.size(); i++) {
            RECLAIM_ADEPT_MEMORY();
            shared_ptr<FilterDAG<real> > f = it->second[i]->copy_sameclass();
            if (initial_check(context, f, false)) {
                string simplify_str = simplify(f);
                if (seen.count(simplify_str) == 0) {
                    seen.insert(simplify_str);
                    pop.push_back(make_shared<GAIndividual<real> >(context, f));
                    if (pop.size() >= params.ga_population) {
                        goto done_init_population;
                    }
                }
            }
        }
    }
    done_init_population:
    if (set_pareto) {
        set_population_pareto(context, pop);
    }
}

template<class real>
vector<double> get_tournament_objective(const vector<shared_ptr<GAIndividual<real> > > &pop_prev, const vector<int> &sub, int tournament_size) {
    if (params.tournament_error) {
        shared_ptr<Pareto<int> > pareto = make_shared<Pareto<int> >();
        for (int i = 0; i < tournament_size; i++) {
            shared_ptr<PARETO_TYPE(real)> current = pop_prev[sub[i]]->pareto;
            for (int j = 0; j < (int) current->L.size(); j++) {
                pareto->add(current->L[j].T, current->L[j].err, i);
            }
        }
        vector<double> ans(tournament_size, 0.0);
        
        for (int i = 0; i < (int) pareto->L.size(); i++) {
            int idx = pareto->L[i].p;
            double err = pareto->L[i].err;
            double err_prev = i > 0 ? pareto->L[i-1].err: err;
            double err_next = (i+1 < (int) pareto->L.size()) ? pareto->L[i+1].err: err;
            err      = MIN(err,      params.pareto_max_error_consider);
            err_prev = MIN(err_prev, params.pareto_max_error_consider);
            err_next = MIN(err_next, params.pareto_max_error_consider);
            
            double err_avg = (fabs(err - err_prev) + fabs(err - err_next)) * 0.5;
            ans[idx] -= err_avg;
        }
        return ans;
    } else {
        vector<double> ans;
        for (int i = 0; i < tournament_size; i++) {
            ans.push_back(pop_prev[sub[i]]->avg_time);
        }
        return ans;
    }
}

template<class real>
int tournament_select(const vector<shared_ptr<GAIndividual<real> > > &pop_prev, int tournament_size=-1) {
    if (params.verbose) {
        printf("tournament_select: Begin\n");
    }
    if (tournament_size < 0) { tournament_size = params.ga_tournament_size; }
    ASSERT(tournament_size < (int) pop_prev.size(), "expected tournament_size < population size");
    vector<int> sub(pop_prev.size(), 0);
    for (int i = 0; i < (int) pop_prev.size(); i++) {
        sub[i] = i;
    }
    
    std::random_shuffle(sub.begin(), sub.end());
    int ibest = -1;
    double objective_best = 1e150;
    if (params.verbose) {
        printf("tournament_size: %d\n", tournament_size);
    }
    vector<double> objectiveL = get_tournament_objective(pop_prev, sub, tournament_size);
    for (int i = 0; i < tournament_size; i++) {
        int ip = sub[i];
        if (params.verbose) {
            printf("tournament_size: %d, objective=%f, avg_time=%f\n", tournament_size, pop_prev[ip]->objective, pop_prev[ip]->avg_time);
        }
        if (objectiveL[i] < objective_best) {
            objective_best = objectiveL[i];
            ibest = ip;
        }
    }
    if (params.verbose) {
        printf("tournament_select: Done\n");
    }
    ASSERT(ibest >= 0, "ibest should be >= 0");
    return ibest;
}

template<class real>
void print_generation(FILTER_CONTEXT context, const vector<shared_ptr<GAIndividual<real> > > &pop, int gen, string msg="") {
    logf("\n-------------------------------------------------------------------------------\n");
    logf("Generation %d%s\n", gen, msg.c_str());
    logf("-------------------------------------------------------------------------------\n");
    for (int i = 0; i < (int) pop.size(); i++) {
        bool ok = false;
        if (pop[i]->pareto) {
            auto pL = pop[i]->pareto->L;
            if (pL.size()) {
                logf("  Mean T: %.2f, Left T: %.2f, E: %.2f, Right T: %.2f, E: %.2f (T_best_max: %.2f), objective: %f, out: %s, simplify: %s\n", pop[i]->avg_time, pL[0].T, pL[0].err, pL[pL.size()-1].T, pL[pL.size()-1].err, context->best_max_time, pop[i]->objective, pop[i]->out_path.c_str(), pop[i]->simplify_str.c_str());
                ok = true;
            }
        }
        if (!ok) {
            logf("Missing Pareto list for %s, %s\n", pop[i]->out_path.c_str(), pop[i]->simplify_str.c_str());
        }
    }
    logf("\n");
}

template<class real>
vector<shared_ptr<GAIndividual<real> > > next_generation(FILTER_CONTEXT context, set<string> &seen, const vector<shared_ptr<GAIndividual<real> > > &pop_prev) {
    vector<shared_ptr<GAIndividual<real> > > ans;
    
    double sum_prob = params.ga_frac_elitism + params.ga_frac_mutate + params.ga_frac_crossover;
    double frac_elitism   = params.ga_frac_elitism   / sum_prob;
    double frac_crossover = params.ga_frac_crossover / sum_prob;
    //double frac_mutate    = params.ga_frac_mutate    / sum_prob;
    
    int n_elitism   = int(params.ga_population * frac_elitism);
    int n_crossover = int(params.ga_population * frac_crossover);
    
    for (int i = 0; i < n_elitism; i++) {
        if (i < (int) pop_prev.size()) {
            ans.push_back(pop_prev[i]);
        }
    }
    
    int crossover_count = 0;
    for (int i = 0; i < n_crossover; i++) {
        shared_ptr<FilterDAG<real> > f;
        for (int j = 0; j < params.ga_crossover_attempts; j++) {
            RECLAIM_ADEPT_MEMORY();
            int ai = tournament_select(pop_prev);
            int bi = tournament_select(pop_prev);
            if (params.verbose) {
                printf("Attempt crossover(%d, %d)\n", ai, bi);
            }
            if (ai == bi) { continue; }
            shared_ptr<FilterDAG<real> > a(pop_prev[ai]->f);
            shared_ptr<FilterDAG<real> > b(pop_prev[bi]->f);
            auto fL = crossover(a, b, seen);
            if (fL.size()) {
                f = fL[0];
                if (params.verbose) {
                    printf("Checking crossover result:\n%s\n\n", f->str().c_str());
                }
                if (initial_check(context, f)) {
                    ans.push_back(make_shared<GAIndividual<real> >(context, f));
                    crossover_count++;
                    break;
                } else {
                    if (params.verbose) {
                        printf("Failed initial check\n");
                    }
                }
            } else {
                if (params.verbose) {
                    printf("Failed to crossover\n");
                }
            }
        }
    }

    int n_mutate = params.ga_population - n_elitism - crossover_count; //int(params.ga_population * frac_mutate);

    for (int i = 0; i < n_mutate; i++) {
        shared_ptr<FilterDAG<real> > f;
        for (int j = 0; j < params.ga_crossover_attempts; j++) {
            RECLAIM_ADEPT_MEMORY();
            f = pop_prev[tournament_select(pop_prev)]->f;
            vector<shared_ptr<FilterDAG<real> > > fL;
            int r = rand() % (params.ga_downsample ? 6: 3);
            if (r == 0) {
                fL = add_mutate(f, seen, false);
            } else if (r == 1) {
                fL = replace_mutate(f, seen);
            } else if (r == 2) {
                fL = delete_mutate(f, seen, false);
            } else if (r >= 3) {
                fL = add_mutate_band(f, seen);
            }
            if (fL.size()) {
                f = fL[0];
                if (initial_check(context, f)) {
                    ans.push_back(make_shared<GAIndividual<real> >(context, f));
                    break;
                } else {
                    if (params.verbose) {
                        printf("Failed initial check\n");
                    }
                }
            } else {
                if (params.verbose) {
                    printf("Failed to mutate\n");
                }
            }
        }
    }
    
    if (int(ans.size()) < params.ga_population) {
        if (params.verbose) {
            printf("next_generation: insufficient population (%d), topping off using combinatoric search\n", int(ans.size()));
        }
        init_population(context, seen, ans, false);            /* Top off using combinatoric search */
    }
    if (params.verbose) {
        printf("next_generation: total %d, elitism: %d, crossover: %d, mutate: %d\n", int(ans.size()), n_elitism, crossover_count, n_mutate);
    }
    
    RECLAIM_ADEPT_MEMORY();
    set_population_pareto(context, ans);
    
    return ans;
}

template<class real>
void main_sample(int argc, char *argv[]) {
    params.flatten = true;
    vector<string> argL = parse_args(argc, argv);
    int n = 100;
    if (argL.size() > 0) { n = std::stoi(argL[0]); }
    int min_size = 3;
    if (argL.size() > 1) { min_size = std::stoi(argL[1]); }
    int max_size = 5;
    if (argL.size() > 2) { max_size = std::stoi(argL[2]); }
    
    params.max_combinatoric = max_size;
    FILTER_MAP(real) filters;
    FILTER_MAP(real) filters_checked;
    
    vector<int> sizes({32, 32});
    Array<real> out(sizes);
    Array<real> in = impulse<real>(sizes[1], sizes[0]);
    Array<double> target(sizes);
    target.clear();
    for (int dx = -2; dx <= 2; dx++) {
        for (int dy = -2; dy <= 2; dy++) {
            target(target.height()/2+dy, target.width()/2+dx) = 1;
        }
    }
    target = normalize_target(target);
    FILTER_CONTEXT context = make_shared<FilterContext<real> >(in, out, target, argc, argv, "");
    
    get_combinatoric_filters(filters, true);

    printf("n=%d, min_size=%d, max_size=%d\n\n", n, min_size, max_size);
    
    printf("---------------------------------------\n");
    printf("Before check\n");
    printf("---------------------------------------\n");
    
    filter_count(filters);
    vector<double> dag_sizes;
    
    for (auto item: filters) {
        //int sz = item.first;
        for (auto filter: item.second) {
            if (initial_check(context, filter, false)) {
                int levels = dag_levels(filter, dag_sizes);
                //if (levels > 0) {
                filters_checked[levels].push_back(filter);
                //}
            }
        }
    }
    
    printf("---------------------------------------\n");
    printf("After check\n");
    printf("---------------------------------------\n");

    filter_count(filters_checked);
    
    PARETO_TYPE(real) pareto;           /* Actually just a list of samples -- not filtered to be the Pareto frontier */
    
    for (auto item: filters_checked) {
        int level = item.first;
        if (level >= 0 && level <= 2) {
            random_shuffle(item.second.begin(), item.second.end());
            int nsub = n/3;
            for (int idx = 0; idx < MIN(nsub, int(item.second.size())); idx++) {
                auto filter = item.second[idx];
                auto fsparse = make_shared<FilterSparse<real> >(filter);
                bool success = false;
                for (int randomize_iter = 0; randomize_iter < 11; randomize_iter++) {
                    RECLAIM_ADEPT_MEMORY();
                    randomize_filter(fsparse, randomize_iter < 10 ? -1: int(fsparse->coeffs().size()), true);
                    if (get_1norm_scale(in, out, fsparse, target) > 1e-8) {
                        success = true;
                        break;
                    }
                }
                //real fscale = 1.0;
                //fit_error(in, out, dynamic_pointer_cast<Filter<real> >(f), target, false, &fscale);
                //fit_error_double_sqrt(in, out, fsparse, target, bool rescale_out=false)
                if (success) {
                    auto point = ParetoPoint<shared_ptr<FilterSparse<real> > >(fsparse->time(), 0.0, fsparse);
                    pareto.L.push_back(point);
                } else {
                    nsub++;
                }
            }
        }
    }
    
    write_pareto(in, out, target, &pareto, false, false, "", true);
    printf("Wrote %d filters to %s\n", int(pareto.L.size()), params.pareto_full.c_str());
}

template<class real>
void main_ga(int argc, char *argv[]) {
    double T0 = wall_time();
    if (argc < 1) {
        usage();
    }
    print_begin_solver(argc, argv);

    string target_filename = string(argv[0]);
    argc--; argv++;
    parse_args(argc, argv);

    if (params.out_arg.size() && params.out_arg[params.out_arg.size()-1] != '/') {
        params.out_arg += '/';
    }
    mkdir(params.out_arg.c_str(), 0777);

    srand(params.seed);
    
//    auto pareto_main = new PARETO_TYPE(real)();
    
    GET_TARGET_IN_OUT(target_filename);
    auto context = make_shared<FilterContext<real> >(in, out, target, argc, argv, target_filename);
    set<string> seen;
    
    vector<shared_ptr<GAIndividual<real> > > pop;              /* Has no Pareto front computed yet */
    
    init_population(context, seen, pop, true);
    
    for (int gen = 1; gen <= params.ga_generations; gen++) {
        if (gen == 1) {
            print_generation(context, pop, 0);
        }
        vector<shared_ptr<GAIndividual<real> > > pop_prev(pop);
        pop = next_generation(context, seen, pop_prev);
        
        print_generation(context, pop, gen);
    }
    
//    write_pareto(in, out, target, pareto_main);
    double T1 = wall_time();
    printf("Total time: %f secs\n", T1-T0);
}

#define APPLY_ONCE(rescale) \
pareto_L[index]->apply({&in}, out, temp); \
if ((rescale) && scaleL.size() == pareto_L.size()) { \
    out = out * real(scaleL[index]); \
}

template<class real>
void main_apply(int argc, char *argv[]) {
    if (argc < 4) {
        usage();
    }
    string in_filename = argv[0];
    string json_filename = argv[1];
    int index = std::stoi(string(argv[2]));
    string out_filename = argv[3];
    
    argc -= 4;
    argv += 4;
    parse_args(argc, argv);
    
    Array<real> in = load_matrix<double>(in_filename);
    ASSERT(in.sizes.size() == 2, "expected 2D target");
    Array<real> out = in;
    
    if (params.target_w > 0) {
        in     = zero_pad_center(in,     params.target_w, params.target_w);
        out    = zero_pad_center(out,    params.target_w, params.target_w);
    }

    vector<double> scaleL;
    vector<shared_ptr<FilterSparse<real> > > pareto_L = json_to_filters<real>(json_filename, &scaleL);
    if (index < 0) {
        index += int(pareto_L.size());
    }
    ASSERT2((unsigned) index < pareto_L.size(), "index out of bounds");

    vector<Array<real> > temp;
    APPLY_ONCE(true);
    
    if (params.time_apply) {
        int iters = 1000;
        APPLY_ONCE(false);
//        array_debug = true;
//        printf("Begin time adouble\n");
        double T0, T1;
        T_fir = 0.0;
        T0 = wall_time();
        for (int i = 0; i < iters; i++) {
            RECLAIM_ADEPT_MEMORY();
//            printf(" ------ Before apply ------\n");
            APPLY_ONCE(false);
//            printf(" ------ After apply ------\n");
        }
        T1 = wall_time();
        printf("Time per iter (adouble): %f\n", (T1-T0)/iters);
        printf("Total time: %f, FIR time: %f\n", T1-T0, T_fir);
        APPLY_ONCE(true);

        Array<double> in_double(in);
        Array<double> out_double(out);
        vector<Array<double> > temp_double;

        pareto_L[index]->apply({&in_double}, out_double, temp_double);

        printf("Begin time double\n");
        T0 = wall_time();
        T_fir = 0.0;
        for (int i = 0; i < iters; i++) {
            RECLAIM_ADEPT_MEMORY();
            pareto_L[index]->apply({&in_double}, out_double, temp_double);
        }
        T1 = wall_time();
        printf("Time per iter (double): %f\n", (T1-T0)/iters);
    }
    
    save_matrix(out, out_filename);
}

template<class real>
void main_solver(int argc, char *argv[]) {
    double T0 = wall_time();

    string solver_cmd = "./solve_kernel solve";
    for (int i = 0; i < argc; i++) {
        solver_cmd += " ";
        solver_cmd += argv[i];
    }

    if (argc < 1) {
        usage();
    }
    print_begin_solver(argc, argv);
    string target_filename = string(argv[0]);
    argc--; argv++;

    FILTER_MAP(real) filters;
    bool init = false;
    if (argc >= 1 && string(argv[0]).size() && string(argv[0])[0] != '-') {
        string json_filename = string(argv[0]);
        argc--; argv++;
        init = true;
        vector<shared_ptr<FilterDAG<real> > > L = json_to_filters_dag<real>(json_filename);
        for (int i = 0; i < (int) L.size(); i++) {
            filters[L[i]->L.size()].push_back(L[i]);
        }
    }
    
    parse_args(argc, argv);

    srand(params.seed);
    
    if (!init) {
        if (params.fixed_topology.size()) {
            get_fixed_topology<real>(filters);
        } else if (params.max_combinatoric > 0) {
            get_combinatoric_filters<real>(filters);
        } else {
            fprintf(stderr, "Use one of: [json.txt] optional arg, -topology, -max_combinatoric\n");
            exit(1);
        }
    }
    
    int nfilter = filter_count<real>(filters);
    int ifilter = 0;
    
    GET_TARGET_IN_OUT(target_filename);

    auto pareto_main = params.resume.size() ? load_pareto(params.resume, in, out, target): make_shared<PARETO_TYPE(real)>();

    for (auto it = filters.begin(); it != filters.end(); it++) {
        for (int i = 0; i < it->second.size(); i++) {
            printf("Filter %d/%d (size %d)\n", ifilter, nfilter, it->first);
            ifilter++;
            
            for (int iter = 0; iter < params.random_iters; iter++) {
                PARETO_TYPE(real) *pareto_sub = pareto_main.get();
                if (params.random_iters > 1) { pareto_sub = new PARETO_TYPE(real)(); }
                shared_ptr<FilterSparse<real> > fsparse;
                if (params.multires) {
                    fsparse = make_shared<FilterSparseMultires<real> >(it->second[i]->copy());
                } else {
                    fsparse = make_shared<FilterSparse<real> >(it->second[i]->copy());
                }
                if (params.verbose >= 2) {
                    printf("fsparse:\n");
                    printf("%s\n\n", fsparse->str().c_str());
                    printf("in:\n");
                    printf("%s\n\n", in.str().c_str());
                    printf("out:\n");
                    printf("%s\n\n", out.str().c_str());
                    printf("target:\n");
                    printf("%s\n\n", target.str().c_str());
                }
                int nlevels = fsparse->nlevels();
                int levels_processed = 0;
                for (int level = nlevels-1; level >= 0; level--) {
                    levels_processed++;
                    if (levels_processed <= params.start_level) {
                        continue;
                    }
                    fsparse->set_level(level);
                    if (params.verbose) {
                        printf("set_level %d (of %d), # coeffs: %d\n", level, nlevels, (int) fsparse->mask.size());
                    }
                    set_pareto_level(pareto_sub, level, in, out, target);
                    if (params.verbose) {
                        printf("# coeffs in Pareto: ");
                        for (int k = 0; k < (int) pareto_sub->L.size(); k++) {
                            printf("%d ", (int) pareto_sub->L[k].p->mask.size());
                        }
                        printf("\n");
                    }
                    optimize_sparse_filter(in, out, fsparse, target, pareto_sub, true, level == nlevels-1);
                    check_pareto(in, out, target, pareto_sub);
                    if (levels_processed >= params.stop_level) {
                        printf("stopping after processing %d levels\n", levels_processed);
                        break;
                    }
                }
                if (params.random_iters > 1) {
                    printf("Joining Pareto\n");
                    pareto_main->join(*pareto_sub);
                    if (params.verbose >= 2) {
                        printf("pareto_sub %d\n", (int) pareto_sub->L.size());
                    }
                    delete pareto_sub;
                }
            }

        }
    }
    
    printf("Writing Pareto\n");
    write_pareto(in, out, target, pareto_main.get(), true, true, solver_cmd);
    double T1 = wall_time();
    printf("Total time: %f secs\n", T1-T0);
}

#endif

