
#include "test_solver.h"
#include <cstdlib>

void test_filters() {
    int w = 7;
    int Kw = 3;
    srand(0);
    
    auto I = Array<double>::random({w, w});
    auto KA = Array<double>::random({Kw, Kw});
    auto KB = Array<double>::random({Kw, Kw});
    
    auto temp = I;
    auto out = I;
    
    auto fA = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(KA));
    auto fB = shared_ptr<FilterFIR<double> >(new FilterFIR<double>(KB));
    auto fA_base = static_pointer_cast<Filter<double> >(fA);
    auto fB_base = static_pointer_cast<Filter<double> >(fB);
    
    out.clear();
    temp.clear();
    
    static vector<Array<double> > temp_images;

    fA->apply({&I}, temp, temp_images);
    fB->apply({&temp}, out, temp_images);
    
    printf("I:\n%s\n\n", I.str().c_str());
    printf("KA:\n%s\n\n", KA.str().c_str());
    printf("KB:\n%s\n\n", KB.str().c_str());
    printf("out:\n%s\n\n", out.str().c_str());
    
    printf("%p\n", (void *) fA->coeffs()[0]);
    
    FilterDAG<double> fL({&fA_base, &fB_base});
    
    for (int iter = 0; iter < 2; iter++) {
        fL.apply({&I}, out, temp_images);
        
        printf("I:\n%s\n\n", I.str().c_str());
        printf("out:\n%s\n\n", out.str().c_str());
    }
}

/* Main test routine */
void test_solver() {
    //test_filters();
    /*
     double T_num  = test_fit<double>();
     double T_auto = test_fit<adouble>();
     printf("Speedup: %f\n", T_num/T_auto);
     */
    
    printf("=======================================================\n");
    printf("Optimize with automatic differentiation:\n");
    printf("=======================================================\n");
    
    test_fit_bfgs<adouble>();
    
//    printf("=======================================================\n");
//    printf("Optimize with numerical differentiation:\n");
//    printf("=======================================================\n");
    
//    test_fit_bfgs<double>();
    
    printf("=======================================================\n");
    printf("Pareto Frontier (Sparse):\n");
    printf("=======================================================\n");
    
    test_fit_sparse<adouble>();
}

void set_out_prefix(const string &a) {
    params.out_arg = a;
    params.pareto      = a + "_pareto_summary.txt";
    params.pareto_full = a + "_pareto_full.txt";
    params.trace       = a + "_trace.json";
}

/* Parse arguments with dashes (e.g. -random 100), and return positional arguments */
vector<string> parse_args(int argc, char *argv[]) {
    vector<string> ans;
    for (int i = 0; i < argc; i++) {
        string s = argv[i];
        if (s.size() >= 1 && s[0] == '-' && (i+1) < argc) {
            string a = argv[i+1];
            if (s == "-random")                           { params.random_iters                   = std::stoi(a); }
            else if (s == "-contract")                    { params.contract_iters                 = std::stoi(a); }
            else if (s == "-expand")                      { params.expand_iters                   = std::stoi(a); }
            else if (s == "-translate")                   { params.translate_iters                = std::stoi(a); }
            else if (s == "-iters")                       { params.translate_iters = params.contract_iters = params.expand_iters = std::stoi(a); }
            else if (s == "-filter_w") {
                params.filter_w               = std::stoi(a);
                if (params.filter_w % 2 != 1) {
                    fprintf(stderr, "expected odd -filter_w\n");
                    exit(1);
                }
            }
            else if (s == "-target_w")                    { params.target_w                       = std::stoi(a); }
            else if (s == "-verbose")                     { params.verbose                        = std::stod(a); }
            else if (s == "-topology")                    { params.fixed_topology                 = a; }
            else if (s == "-pareto")                      { params.pareto                         = a; }
            else if (s == "-pareto_full")                 { params.pareto_full                    = a; }
            else if (s == "-truncate_row_col")            { params.truncate_row_col               = std::stoi(a); }
            else if (s == "-multires")                    { params.multires                       = std::stoi(a); }
            else if (s == "-stop_level")                  { params.stop_level                     = std::stoi(a); }
            else if (s == "-start_level")                 { params.start_level                    = std::stoi(a); }
            else if (s == "-feature_coeff")               { params.feature_coeff                  = a; }
            else if (s == "-solve_initial_maxiters")      { params.solve_initial_maxiters         = std::stoi(a); }
            else if (s == "-solve_initial_epsilon")       { params.solve_initial_epsilon          = std::stod(a); }
            else if (s == "-solve_add_maxiters")          { params.solve_add_maxiters             = std::stoi(a); }
            else if (s == "-solve_add_epsilon")           { params.solve_add_epsilon              = std::stod(a); }
            else if (s == "-reoptim_add")                 { params.reoptim_add                    = std::stoi(a); }
            else if (s == "-reoptim_thresh")              { params.reoptim_thresh                 = std::stod(a); }
            else if (s == "-contract_is_random")          { params.contract_is_random             = std::stoi(a); }
            else if (s == "-expand_is_random")            { params.expand_is_random               = std::stoi(a); }
            else if (s == "-expand_recheck")              { params.expand_recheck                 = std::stoi(a); }
            else if (s == "-contract_recheck")            { params.contract_recheck               = std::stoi(a); }
            else if (s == "-seed") {
                params.seed = std::stoi(a);
                srand(params.seed);
                std::srand(params.seed);
            } else if (s == "-translate_step")              { params.translate_step                 = std::stod(a); }
            else if (s == "-translate_all")               { params.translate_all                  = bool(std::stoi(a)); }
            else if (s == "-translate_some")              { params.translate_some                 = bool(std::stoi(a)); }
            else if (s == "-translate_duplicate")         { params.translate_duplicate            = bool(std::stoi(a)); }
            else if (s == "-translate_duplicate_outside") { params.translate_duplicate_outside    = bool(std::stoi(a)); }
            else if (s == "-translate_cross")             { params.translate_cross                = bool(std::stoi(a)); }
            else if (s == "-translate_border")            { params.translate_border               = bool(std::stoi(a)); }
            else if (s == "-translate_copy")              { params.translate_copy                 = bool(std::stoi(a)); }
            else if (s == "-random_order")                { params.random_order                   = bool(std::stoi(a)); }
            else if (s == "-visit_one")                   { params.visit_one                      = bool(std::stoi(a)); }
            else if (s == "-scale_correct")               { params.scale_correct                  = bool(std::stoi(a)); }
            else if (s == "-search_shifts")               { params.search_shifts                  = bool(std::stoi(a)); }
            else if (s == "-shift_xmin")                  { params.shift_xmin                     = std::stoi(a); }
            else if (s == "-shift_ymin")                  { params.shift_ymin                     = std::stoi(a); }
            else if (s == "-shift_xmax")                  { params.shift_xmax                     = std::stoi(a); }
            else if (s == "-shift_ymax")                  { params.shift_ymax                     = std::stoi(a); }
            else if (s == "-check_often")                 { params.check_often                    = bool(std::stoi(a)); }
            else if (s == "-trace")                       { params.trace                          = a; }
            else if (s == "-uniform_error")               { params.uniform_error                  = bool(std::stoi(a)); }
            else if (s == "-uniform_error")               { params.uniform_error                  = bool(std::stoi(a)); }
            else if (s == "-recalc_solve")                { params.recalc_solve                   = bool(std::stoi(a)); }
            else if (s == "-symmetry")                    { params.symmetry                       = bool(std::stoi(a)); }
            else if (s == "-symmetry_w_min")              { params.symmetry_w_min                 = std::stod(a); }
            else if (s == "-symmetry_w_max")              { params.symmetry_w_max                 = std::stod(a); }
            else if (s == "-random_init")                 { params.random_init                    = bool(std::stoi(a)); }
            else if (s == "-random_count")                { params.random_count                   = std::stoi(a); }
            else if (s == "-max_of_calls")                { params.max_of_calls                   = std::stoi(a); }
            else if (s == "-adaptive")                    { params.adaptive                       = bool(std::stoi(a)); }
            else if (s == "-adaptive_iters")              { params.adaptive_iters                 = std::stoi(a); }
            else if (s == "-try_separable")               { params.try_separable                  = bool(std::stoi(a)); }
            else if (s == "-merge_base")                  { params.merge_base                     = bool(std::stoi(a)); }
            else if (s == "-multisize")                   { params.multisize                      = bool(std::stoi(a)); }
            else if (s == "-multisize_sizes")             { params.multisize_sizes                = std::stoi(a); }
            else if (s == "-quantize")                    { params.quantize                       = bool(std::stoi(a)); }
            else if (s == "-quantize_count")              { params.quantize_count                 = std::stoi(a); }
            else if (s == "-quantize_until")              { params.quantize_until                 = std::stoi(a); }
            else if (s == "-terminate_early")             { params.terminate_early                = bool(std::stoi(a)); }
            else if (s == "-terminate_epsilon")           { params.terminate_epsilon              = std::stod(a); }
            else if (s == "-terminate_iters")             { params.terminate_iters                = std::stoi(a); }
            else if (s == "-mask_samples")                { params.mask_samples                   = std::stoi(a); }

            else if (s == "-max_combinatoric")            { params.max_combinatoric               = std::stoi(a); }
            else if (s == "-ga_downsample")               { params.ga_downsample                  = bool(std::stoi(a)); }
            else if (s == "-ga_always_downsample")        { params.ga_always_downsample           = bool(std::stoi(a)); }
            else if (s == "-ga_generations")              { params.ga_generations                 = std::stoi(a); }
            else if (s == "-ga_population")               { params.ga_population                  = std::stoi(a); }
            else if (s == "-ga_frac_elitism")             { params.ga_frac_elitism                = std::stod(a); }
            else if (s == "-ga_frac_crossover")           { params.ga_frac_crossover              = std::stod(a); }
            else if (s == "-ga_frac_mutate")              { params.ga_frac_mutate                 = std::stod(a); }
            else if (s == "-ga_tournament_size")          { params.ga_tournament_size             = std::stoi(a); }

            else if (s == "-pareto_error_thresh")         { params.pareto_error_thresh            = std::stod(a); }
            else if (s == "-pareto_min_points")           { params.pareto_min_points              = std::stoi(a); }
            else if (s == "-pareto_min_points_error")     { params.pareto_min_points_error        = std::stod(a); }

            else if (s == "-antialias")                   { params.antialias                      = bool(std::stoi(a)); }
            else if (s == "-antialias_subsample")         { params.antialias_subsample            = std::stoi(a); }
            else if (s == "-antialias_levels")            { params.antialias_levels               = std::stoi(a); }
            else if (s == "-antialias_interval")          { params.antialias_interval             = bool(std::stoi(a)); }

            else if (s == "-out")                         { set_out_prefix(a); }
            else if (s == "-in")                          { params.in                             = a; }

            else if (s == "-preconvolve")                 { params.preconvolve = a; }

            else if (s == "-prefilter_sigma")             { params.prefilter_sigma        = std::stod(a); }
            else if (s == "-prefilter_size")              { params.prefilter_size         = std::stoi(a); }
            
            else if (s == "-time_apply")                  { params.time_apply                     = bool(std::stoi(a)); }
            else if (s == "-flatten")                     { params.flatten                        = bool(std::stoi(a)); }
            else if (s == "-retain_all")                  { params.retain_all                     = bool(std::stoi(a)); }

            else if (s == "-feature")                     { params.feature                        = std::stoi(a); }
            else if (s == "-resume")                      { params.resume                         = a; }
            else if (s == "-weights")                     { params.weights                        = a; }
            else if (s == "-ignore_boundary")             { params.ignore_boundary                = std::stoi(a); }
            else if (s == "-scale_1norm")                 { params.scale_1norm                    = bool(std::stoi(a)); }

            else if (s == "-max_time")                    { params.max_time                          = std::stod(a); }
            else if (s == "-vh_mode")                     { params.vh_mode                           = bool(std::stoi(a)); }
            else if (s == "-seed_convpyr")                { params.seed_convpyr                      = bool(std::stoi(a)); }

            else { fprintf(stderr, "Error parsing switch %s\n", s.c_str()); usage(); }
            i++;
        } else if (s.size() >= 1 && s[0] != '-') {
            ans.push_back(s);
        } else {
            fprintf(stderr, "Error parsing argument %s\n", s.c_str());
            usage();
        }
    }
    return ans;
}

void print_begin_solver(int argc, char *argv[]) {
    printf("========================================================\n");
    printf("Main solver\n");
    printf("========================================================\n\n");
    printf("%s", "% ./solve_kernel solve ");
    for (int i = 0; i < (int) argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    string git_rev = exec_output("git rev-parse --short HEAD");
    printf("Git revision %s\n", git_rev.c_str());
}

void run_commands(const vector<string> &L) {
/*    for (int i = 0; i < (int) L.size(); i++) {
        printf("%s\n", L[i].c_str());
        system(L[i].c_str());
    }*/
    if (!L.size()) { return; }
    string filename = params.out_arg + "run.sh";
    FILE *f = fopen(filename.c_str(), "wt");
    for (int i = 0; i < (int) L.size(); i++) {
        fprintf(f, "%s\n", L[i].c_str());
    }
    fclose(f);
    char buf[1024];
    sprintf(buf, "python parallel.py %s", filename.c_str());
    system(buf);
}
