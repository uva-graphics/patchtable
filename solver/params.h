
#ifndef _params_h
#define _params_h

#include <string>
#include <vector>
#include "array.h"
using std::string;
using std::vector;

#define DEFAULT_SPARSEOPTIM_ITERS 2000

#define DEFAULT_SOLVER_MAXITERS 0
#define DEFAULT_SOLVER_EPSILON  1e-5

/* Parameters to solver */
class Params {
    vector<double> feature_coeff_list_internal;
    Array<double> weights_array_internal;
    
public:
    bool truncate_taylor_approx;
    int optimize_top_k;
    int random_iters;
    int contract_iters;
    int translate_iters;
    int expand_iters;
    bool init_correct_sparsity;
    bool expand_is_random;
    bool contract_is_random;
    int kernel_function;
    int filter_w;
    int target_w;                           /* If positive, crops or pads the target to size target_w x target_w */
    int max_combinatoric;
    bool truncate_row_col;
    string fixed_topology;
    string out_arg;
    string pareto;
    string pareto_full;
    bool multires;
    
    double verbose;
    int start_level;
    int stop_level;
    
    int solve_initial_maxiters;
    double solve_initial_epsilon;
    int solve_add_maxiters;
    double solve_add_epsilon;
    bool reoptim_add;                       /* Re-optimize when adding? */
    double reoptim_thresh;                  /* If distance is <= thresh then reoptimize and add */
    
    bool contract_recheck;
    bool expand_recheck;
    bool translate_all;
    bool translate_some;
    bool translate_duplicate;
    bool translate_duplicate_outside;
    bool translate_border;
    bool translate_cross;
    bool translate_copy;
    double translate_step;
    
    int seed;
    bool random_order;
    bool visit_one;
    bool scale_correct;
    
    bool search_shifts;
    int shift_xmin, shift_xmax, shift_ymin, shift_ymax;
    
    bool check_often;
    string trace;
    bool uniform_error;
    int uniform_error_samples;
    bool recalc_solve;
//    bool spatial_weights;
//    double spatial_sigma_frac;
    
    bool symmetry;                          /* Enforce symmetry? */
    bool symmetry_h;                        /* Detected H symmetry */
    bool symmetry_v;                        /* Detected V symmetry */
    double symmetry_w_min;                  /* Symmetry term is multiplied by w_min * pow(w_max/w_min, t), t=iter/max_iters */
    double symmetry_w_max;
    bool allow_symmetry;
    int symmetry_reset_iters;
    
    bool random_init;
    int random_count;
    
    int max_of_calls;
    bool adaptive;
    int adaptive_iters;
    
    bool try_separable;                     /* Try all horizontal, all vertical, unit impulse, for FIR filters */
    bool merge_base;                        /* Merge with base (single FIR) pareto before writing */
    
    bool multisize;                         /* Solve using multiple sizes, starting with small sizes. */
    int multisize_sizes;                    /* Number of sizes */
    
    bool quantize;                          /* Quantize error */
    int quantize_count;
    int quantize_until;
    
    bool terminate_early;
    double terminate_epsilon;
    int terminate_iters;
    
    int mask_samples;                       /* If positive, sample each mask at most n times */
    int max_resample;
    
    double max_time;                        /* Max time in minutes for solver */
    
    bool ga_downsample;                     /* Permit downsample in GA */
    bool ga_always_downsample;              /* Always downsample in GA */
    int ga_crossover_attempts;
    int ga_generations;
    int ga_population;
    
    double ga_frac_elitism;
    double ga_frac_crossover;
    double ga_frac_mutate;
    int ga_tournament_size;
    
    double pareto_error_thresh;             /* avg_time in the G.A. is sorted by area under the curve from [0, pareto_error_thresh]  */
    double pareto_min_points_error;
    int pareto_min_points;                  /* minimum number of points with error <= pareto_min_points_error */
    double pareto_max_error_consider;
    
    bool antialias;                         /* Whether to antialias */
    int antialias_subsample;                /* If positive, take n samples from total set of antialiasing offsets. */
    int antialias_levels;
    bool antialias_interval;
    
    string in;                              /* Load input image (2D matrix) */
    string preconvolve;
    
    double prefilter_sigma;
    int prefilter_size;
    
    bool time_apply;
    
    bool flatten;
    bool retain_all;
    int feature;
    
    string resume;
    
    string feature_coeff;
    int ignore_boundary;
    
    double target_orig_norm;
    bool scale_1norm;
    string weights;                         /* 2D Matrix of weights (same size as target) */
    
    double time_limit;
    bool vh_mode;
    bool seed_convpyr;
    bool tournament_error;
    
    vector<double> *feature_coeff_list();
    Array<double> *weights_array();
    
    Params()
    :truncate_taylor_approx(false),         /* Use (linear) Taylor approx when truncating / increasing sparsity */
    optimize_top_k(1000*1000),
    random_iters(1),
    contract_iters(DEFAULT_SPARSEOPTIM_ITERS),
    translate_iters(DEFAULT_SPARSEOPTIM_ITERS),
    expand_iters(DEFAULT_SPARSEOPTIM_ITERS),
    init_correct_sparsity(false),
    expand_is_random(true),
    contract_is_random(true),
    kernel_function(0),
    filter_w(7),
    target_w(-1),
    max_combinatoric(4),
    truncate_row_col(false),
    fixed_topology(""),
    out_arg("output"),
    pareto("pareto.txt"),
    pareto_full("pareto_full.txt"),
    multires(false),
    verbose(0),                             /* Verbose level: 0, 1, 2, 3. */
    start_level(0),
    stop_level(1000),
    solve_initial_maxiters(5),    // DEFAULT_SOLVER_MAXITERS
    solve_initial_epsilon(1e-12),
    solve_add_maxiters(10),
    solve_add_epsilon(1e-12),
    reoptim_add(true),
    reoptim_thresh(1.0),
    contract_recheck(true),
    expand_recheck(false),
    translate_all(false),
    translate_some(false),
    translate_duplicate(true),
    translate_duplicate_outside(false),
    translate_border(true),       // false
    translate_cross(false),
    translate_copy(false),
    translate_step(1),
    seed(0),
    random_order(true),
    visit_one(false),
    scale_correct(true),         // TODO: This should probably be true
    search_shifts(false),
    shift_xmin(-1),
    shift_xmax(1),
    shift_ymin(-1),
    shift_ymax(1),
    check_often(false),
    trace("trace.json"),
    uniform_error(false),
    uniform_error_samples(100),
    recalc_solve(false),
    symmetry(false),
    symmetry_h(false),
    symmetry_v(false),
    symmetry_w_min(0.01),
    symmetry_w_max(100.0),
    allow_symmetry(true),        // Internal use only
    symmetry_reset_iters(100),
    random_init(true),
    random_count(500),
    max_of_calls(-1),
    adaptive(false),
    adaptive_iters(100),          // Freeze a Pareto point after this many iters
    try_separable(false),
    merge_base(false),
    multisize(false),
    multisize_sizes(5),
    quantize(true),
    quantize_count(1000),
    quantize_until(-1),
    terminate_early(true),
    terminate_epsilon(1e-3),
    terminate_iters(150),
    mask_samples(1),
    max_resample(2),
    max_time(-1),
    
    ga_downsample(true),
    ga_always_downsample(false),
    ga_crossover_attempts(10),
    ga_generations(10),
    ga_population(30),
    
    ga_frac_elitism(0.33),
    ga_frac_crossover(0.33),
    ga_frac_mutate(0.33),
    ga_tournament_size(3),
    
    pareto_error_thresh(0.3),
    pareto_min_points_error(0.15),
    pareto_min_points(1),
    pareto_max_error_consider(0.5),
    
    antialias(false),
    antialias_subsample(-1),
    antialias_levels(0),
    antialias_interval(false),
    
    in(""),
    preconvolve(""),
    
    prefilter_sigma(1.335), //(1.2),
    prefilter_size(7), //(5)
    time_apply(false),
    flatten(false),
    
    retain_all(false),
    feature(1),
    resume(""),
    
    feature_coeff("feature_coeff.txt"),
    
    ignore_boundary(0),
    
    target_orig_norm(1.0),
    scale_1norm(false),
    weights(""),
    
    vh_mode(false),
    seed_convpyr(false),
    tournament_error(true)
//    spatial_weights(true),
//    spatial_sigma_frac(1.0/3.0)
    {}
};

extern Params params;

#endif
