
#include "patchtable.h"

/* ------------------------------------------------------------------------------------
   Parameters, which can be initialized from command-line switches
   ------------------------------------------------------------------------------------ */

#define DEFAULT_PARTITION_STEP 1

map<string, string> parse_switches(int argc, const char **argv) {
    map<string, string> ans;
    for (int i = 1; i < argc; i++) {
        if (strlen(argv[i]) && argv[i][0] == '-' && i+1 < argc) {
            ans[argv[i]] = argv[i+1];
            i += 1;
        }
    }
 //   printf("switches size: %d\n", ans.size());
    return ans;
}

map<string, string> parse_switches(int argc, char **argv) {
    return parse_switches(argc, (const char **) argv);
}

void PatchTableParams::set_defaults() {
    randomize_dt = false;
    parallel_dt = false;
    limit = 0.1e6;          // Was 4e6;
    threads = 1;
    dt_threads = 1;
    convert_colorspace = true;
    colorspace = PATCHTABLE_COLORSPACE_YUV;
    regular_grid = false;
    verbose = 1;
    table_knn = 1;
    populate_nearest = true;
    partition_step = DEFAULT_PARTITION_STEP;    // Was 1
    ndims = 6;
    nchroma = 1;
    ntables = 1;
    populate_exptime = false;
    init_random = false;
    dt_mode = DT_MODE_MANHATTAN;
    dt_algo = DT_ALGO_PROP;
    dt_iters = -1;                  /* Unlimited iters */
    query_step = 1;
    run_dt = true;         // Was true
    calc_exact_dist = false;
    do_prop = true;
    prop_y_step = 1;
    do_rs = true;
    allowed_index = 1;
    spatial = 0;
    prob_rs = 1.0;
    do_calc_dist = true;
    patch_w = 8;
    kcoherence_step = 1;
    coherence_spatial = 0.0;
    coherence_temporal = 0.0;
    table_ratio = 2.0;
    treecann_agrid = 1;
    treecann_bgrid = 1;
    recalc_dist_temporal = true;
    is_descriptor = false;
    check_existing = false;
    descriptor_padded = false;
    lookup_algo = LOOKUP_ALGO_TABLE;
    flann_build_step = 1;
    flann_trees = 4;
    flann_checks = 1;
    flann_eps = 500000;
    treecann_eps = 0;
    flann_reset_step = 10;
    query_iters = 2;
    dt_knn = 1;
    prop_iters = 0;
    enrich_limit = 0.1e6;
    kcoherence = 0;
    kcoherence_min_dist = TABLE_DEFAULT_MIN_DIST;
    pm_iters = 5;
    pm_knn = 0;
    pm_min_dist = TABLE_DEFAULT_MIN_DIST;
    pm_enrich = false;
    kcoherence_algo = KCOHERENCE_ALGO_FLANN;
    pca_samples = 1000;
    dim_algo = TABLE_DIM_ALGO_WH;
    pca_mean = false;
    pca_subtract_mean = true;
    auto_filter = false;
    kdtree_add_unique = false;
    populate_random = true;
    kcoherence_iter = -1;
    kcoherence_enrich = 0;
    grid_ndims = -1;
    prop_dist = 0;
    incoherence = false;
    incoherence_min_dist = TABLE_DEFAULT_MIN_DIST;
    save_kcoherence = false;
    triangle_factor = 2.0;
    sanitize_input = true;
    do_table_lookup = true;
    init_random = false;

#if TABLE_CLUSTER_KMEANS
    cluster_kmeans = false;
    cluster_count = 1024;
    cluster_limit = limit;
#endif

#if TABLE_PRODUCT_QUANTIZE
    product_quantize = false;
    product_quantize_log = false;
    product_quantize_dims = 2;
    product_quantize_mapn = 128;
    product_quantize_knn = 4;
    product_quantize_dt_all = false;
#endif

#if (TABLE_CLUSTER_KMEANS||TABLE_PRODUCT_QUANTIZE)
    kmeans_attempts = 1;
    kmeans_max_iters = 10;
    kmeans_eps = 1.0;
#endif
    
#if TABLE_ENABLE_ANN
    ann = false;
    ann_agrid = 4;
    ann_bgrid = 4;
    ann_eps = 3.0;
#endif
#if TABLE_RS_TABLE
    rs_table = 0;
#endif
    set_speed(3);
}

void PatchTableParams::set_speed(int i) {
    if      (i == 0) { set_from_switches(string("-ndims 8 -nchroma 1 -limit 10.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 1 -kcoherence 20 -kcoherence_iter -1 -prop_iters 2 -spatial 1")); }
    else if (i == 1) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 1 -kcoherence 10 -kcoherence_iter 0 -prop_iters 2 -spatial 1")); }
    else if (i == 2) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 1 -kcoherence 10 -kcoherence_iter 0 -prop_iters 2 -spatial 0")); }
    else if (i == 3) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 0 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 1 -kcoherence 2 -kcoherence_iter -1 -prop_iters 2 -spatial 0")); }
    else if (i == 4) { set_from_switches(string("-ndims 7 -nchroma 1 -limit 10.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 2 -kcoherence 20 -kcoherence_iter 0 -prop_iters 1 -spatial 1")); }
    else if (i == 5) { set_from_switches(string("-ndims 8 -nchroma 1 -limit 100.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 2 -kcoherence 4 -kcoherence_iter 0 -prop_iters 1 -spatial 1")); }
    else if (i == 6) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 0 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 2 -kcoherence 2 -kcoherence_iter -1 -prop_iters 1 -spatial 0")); }
    else if (i == 7) { set_from_switches(string("-ndims 7 -nchroma 1 -limit 1.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 4 -kcoherence 3 -kcoherence_iter 0 -prop_iters 1 -spatial 1")); }
    else if (i == 8) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 6 -kcoherence 5 -kcoherence_iter 0 -prop_iters 1 -spatial 0")); }
    else if (i == 9) { set_from_switches(string("-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 0 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 10 -kcoherence 5 -kcoherence_iter 0 -prop_iters 1 -spatial 0")); }
    else { fprintf(stderr, "Error: set_speed(i), i=%d, expected integer in 0...9\n"); exit(1); }
}

void PatchTableParams::set_from_switches(string s) {
    vector<string> L = str_split(s, ' ');
    vector<const char *> argv;
    for (int i = 0; i < (int) L.size(); i++) {
        if (L[i].size() && !(L[i].size() == 1 && L[i][0] == ' ')) {
            argv.push_back(L[i].c_str());
        }
    }
    set_from_switches(argv.size(), &argv[0]);
}

void PatchTableParams::set_from_switches(int argc, char **argv) {
    set_from_switches(argc, (const char **) argv);
}

void PatchTableParams::set_from_switches(int argc, const char **argv) {
    map<string, string> switches = parse_switches(argc, argv);
    if (switches.count("-speed"))  { set_speed(atoi(switches["-speed"].c_str())); }

    if (switches.count("-limit")) { limit = int(1000*1000*atof(switches["-limit"].c_str())); }
    if (switches.count("-enrich_limit")) { enrich_limit = int(1000*1000*atof(switches["-enrich_limit"].c_str())); }
    
    if (switches.count("-regular_grid")) { regular_grid = atoi(switches["-regular_grid"].c_str()); }
#if TABLE_OPTIMIZE_REGULAR
    regular_grid = true;
#endif

    if (switches.count("-table_knn")) { table_knn = atoi(switches["-table_knn"].c_str()); }
    if (switches.count("-randomize_dt")) { randomize_dt = atoi(switches["-randomize_dt"].c_str()); }
    if (switches.count("-parallel_dt")) { parallel_dt = atoi(switches["-parallel_dt"].c_str()); }
    if (switches.count("-dt_threads")) { dt_threads = atoi(switches["-dt_threads"].c_str()); }
#if TABLE_OPENMP
    if (switches.count("-threads"))  { threads = atoi(switches["-threads"].c_str()); }
#endif
    if (switches.count("-verbose")) { verbose = atoi(switches["-verbose"].c_str()); }
    if (switches.count("-populate_nearest")) { populate_nearest = atoi(switches["-populate_nearest"].c_str()); }
    if (switches.count("-populate_random")) { populate_random = bool(atoi(switches["-populate_random"].c_str())); }
    if (switches.count("-partition_step")) { partition_step = atoi(switches["-partition_step"].c_str()); }

    if (switches.count("-ndims")) { ndims = atoi(switches["-ndims"].c_str()); }
    if (switches.count("-nchroma")) { nchroma = atoi(switches["-nchroma"].c_str()); }
    if (switches.count("-ntables")) { ntables = atoi(switches["-ntables"].c_str()); }

    if (switches.count("-nslices")) {
        vector<string> nL = str_split(switches["-nslices"], ',');
        ndims = nL.size();
        nslices.resize(ndims);
        ASSERT2(nL.size() == nslices.size(), "expected -nslices size to match number of dims");
        for (int i = 0; i < (int) nL.size(); i++) {
            nslices[i] = atoi(nL[i].c_str());
        }
        printf("Parsed nslices: ");
        for (int i = 0; i < (int) nL.size(); i++) {
            printf("%d ", nslices[i]);
        }
        printf("\n");
    }
    
    if (switches.count("-run_dt")) { run_dt = bool(atoi(switches["-run_dt"].c_str())); }
    if (switches.count("-do_prop")) { do_prop = atoi(switches["-do_prop"].c_str()); }
    if (switches.count("-prop_y_step")) { prop_y_step = atof(switches["-prop_y_step"].c_str()); }
    if (switches.count("-do_rs"))   { do_rs = atoi(switches["-do_rs"].c_str()); }
    if (switches.count("-spatial"))  { spatial = atoi(switches["-spatial"].c_str()); }
    if (switches.count("-prob_rs"))  { prob_rs = atof(switches["-prob_rs"].c_str()); }
#if TABLE_RS_TABLE
    if (switches.count("-rs_table"))  { rs_table = atoi(switches["-rs_table"].c_str()); printf("set rs_table=%d\n", rs_table); }
#endif

    if (switches.count("-load"))  { load_filename = switches["-load"]; }
    if (switches.count("-save"))  { save_filename = switches["-save"]; }
    if (switches.count("-check_existing")) { check_existing = bool(atoi(switches["-check_existing"].c_str())); }
        
    if (switches.count("-patch_w"))  { patch_w = atoi(switches["-patch_w"].c_str()); }

    if (switches.count("-pca_samples"))  { pca_samples = atoi(switches["-pca_samples"].c_str()); }
    if (switches.count("-pca_mean"))  { pca_mean = bool(atoi(switches["-pca_mean"].c_str())); }
    if (switches.count("-pca_subtract_mean"))  { pca_subtract_mean = bool(atoi(switches["-pca_subtract_mean"].c_str())); }

#if TABLE_ENABLE_ANN
    if (switches.count("-ann"))  { ann = bool(atoi(switches["-ann"].c_str())); }
    if (switches.count("-ann_agrid"))  { ann_agrid = atoi(switches["-ann_agrid"].c_str()); }
    if (switches.count("-ann_bgrid"))  { ann_bgrid = atoi(switches["-ann_bgrid"].c_str()); }
    if (switches.count("-ann_eps"))  { ann_eps = atof(switches["-ann_eps"].c_str()); }
#endif

#if TABLE_CLUSTER_KMEANS
    if (switches.count("-cluster_kmeans"))  { cluster_kmeans = bool(atoi(switches["-cluster_kmeans"].c_str())); }
    if (switches.count("-cluster_count"))  { cluster_count = atoi(switches["-cluster_count"].c_str()); }
    if (switches.count("-cluster_limit")) { cluster_limit = int(1000*1000*atof(switches["-cluster_limit"].c_str())); }
#endif

    if (switches.count("-dt_knn"))  { dt_knn = atoi(switches["-dt_knn"].c_str()); }

    if (switches.count("-prop_iters"))  { prop_iters = atoi(switches["-prop_iters"].c_str()); }
    
    if (switches.count("-dt_iters"))  { dt_iters = atoi(switches["-dt_iters"].c_str()); }
    if (switches.count("-query_step"))  { query_step = atoi(switches["-query_step"].c_str()); }
    if (switches.count("-coherence_spatial"))  { coherence_spatial = atof(switches["-coherence_spatial"].c_str()); }
    if (switches.count("-coherence_temporal"))  { coherence_temporal = atof(switches["-coherence_temporal"].c_str()); }
    if (switches.count("-recalc_dist_temporal"))  { recalc_dist_temporal = bool(atoi(switches["-recalc_dist_temporal"].c_str())); }

    if (switches.count("-dt_algo")) {
        string mode = switches["-dt_algo"];
        if (mode == string("raster")) { dt_algo = DT_ALGO_RASTER; }
        else if (mode == string("prop")) { dt_algo = DT_ALGO_PROP; }
        else if (mode == string("brute")) { dt_algo = DT_ALGO_BRUTE; }
        else if (mode == string("kdtree")) { dt_algo = DT_ALGO_KDTREE; }
        else if (mode == string("downsample")) { dt_algo = DT_ALGO_DOWNSAMPLE; }
        else if (mode == string("hybrid")) { dt_algo = DT_ALGO_HYBRID; }
		else if (mode == string("euclidean")) { dt_algo = DT_ALGO_EUCLIDEAN; }
        else { fprintf(stderr, "invalid dt_algo\n"); exit(1); }
    }

    if (switches.count("-colorspace")) {
        string mode = switches["-colorspace"];
        if (mode == string("yuv")) { colorspace = PATCHTABLE_COLORSPACE_YUV; }
        else if (mode == string("lab")) { colorspace = PATCHTABLE_COLORSPACE_LAB; }
        else if (mode == string("rgb")) { convert_colorspace = false; }
        else { fprintf(stderr, "unsupported colorspace\n"); exit(1); }
    }

    if (switches.count("-lookup_algo")) {
        string mode = switches["-lookup_algo"];
        if (mode == string("table")) { lookup_algo = LOOKUP_ALGO_TABLE; }
        else if (mode == string("brute")) { lookup_algo = LOOKUP_ALGO_BRUTE; }
        else if (mode == string("kdtree")) { lookup_algo = LOOKUP_ALGO_KDTREE; }
        else if (mode == string("pm")) { lookup_algo = LOOKUP_ALGO_PM; }
        else if (mode == string("treecann")) { lookup_algo = LOOKUP_ALGO_TREECANN; }
        else { fprintf(stderr, "unsupported lookup_algo\n"); exit(1); }
    }

    if (switches.count("-kcoherence_algo")) {
        string mode = switches["-kcoherence_algo"];
        if (mode == string("flann")) { kcoherence_algo = KCOHERENCE_ALGO_FLANN; }
        else if (mode == string("pm")) { kcoherence_algo = KCOHERENCE_ALGO_PM; }
        else if (mode == string("ann")) { kcoherence_algo = KCOHERENCE_ALGO_ANN; }
        else if (mode == string("treecann")) { kcoherence_algo = KCOHERENCE_ALGO_TREECANN; }
        else { fprintf(stderr, "unsupported lookup_algo\n"); exit(1); }
    }

    if (switches.count("-dim_algo")) {
        string mode = switches["-dim_algo"];
        if (mode == string("pca")) { dim_algo = TABLE_DIM_ALGO_PCA; }
        else if (mode == string("wh")) { dim_algo = TABLE_DIM_ALGO_WH; }
        else { fprintf(stderr, "unsupported dim_algo\n"); exit(1); }
    }
    
    if (switches.count("-flann_trees"))  { flann_trees = atoi(switches["-flann_trees"].c_str()); }
    if (switches.count("-flann_checks"))  { flann_checks = atoi(switches["-flann_checks"].c_str()); }
    if (switches.count("-flann_build_step"))  { flann_build_step = atoi(switches["-flann_build_step"].c_str()); }
    if (switches.count("-flann_eps"))  { flann_eps = atof(switches["-flann_eps"].c_str()); }
    if (switches.count("-treecann_eps"))  { treecann_eps = atof(switches["-treecann_eps"].c_str()); }
    if (switches.count("-flann_reset_step"))  { flann_reset_step = atoi(switches["-flann_reset_step"].c_str()); }
    if (switches.count("-query_iters"))  { query_iters = atoi(switches["-query_iters"].c_str()); }
    if (switches.count("-compare"))  { compare = switches["-compare"].c_str(); }

    if (switches.count("-treecann_agrid"))  { treecann_agrid = atoi(switches["-treecann_agrid"].c_str()); }
    if (switches.count("-treecann_bgrid"))  { treecann_bgrid = atoi(switches["-treecann_bgrid"].c_str()); }

    if (switches.count("-kcoherence"))  { kcoherence = atoi(switches["-kcoherence"].c_str()); }
    if (switches.count("-kcoherence_iter"))  { kcoherence_iter = atoi(switches["-kcoherence_iter"].c_str()); }
    if (switches.count("-kcoherence_min_dist"))  { kcoherence_min_dist = atoi(switches["-kcoherence_min_dist"].c_str()); }
    if (switches.count("-kcoherence_step"))  { kcoherence_step = atoi(switches["-kcoherence_step"].c_str()); }
    if (switches.count("-kcoherence_enrich"))  { kcoherence_enrich = atoi(switches["-kcoherence_enrich"].c_str()); }
    if (switches.count("-incoherence"))  { incoherence = bool(atoi(switches["-incoherence"].c_str())); }
    if (switches.count("-incoherence_min_dist"))  { incoherence_min_dist = atoi(switches["-incoherence_min_dist"].c_str()); }
    if (switches.count("-save_kcoherence"))  { save_kcoherence = bool(atoi(switches["-save_kcoherence"].c_str())); }
    if (switches.count("-triangle_factor"))  { triangle_factor = atof(switches["-triangle_factor"].c_str()); }

    if (switches.count("-init_random"))  { init_random = bool(atoi(switches["-init_random"].c_str())); }

    if (switches.count("-pm_iters"))  { pm_iters = atoi(switches["-pm_iters"].c_str()); }
    if (switches.count("-pm_enrich"))  { pm_enrich = bool(atoi(switches["-pm_enrich"].c_str())); }

    if (switches.count("-sanitize_input"))  { sanitize_input = bool(atoi(switches["-sanitize_input"].c_str())); }
    if (switches.count("-do_table_lookup"))  { do_table_lookup = bool(atoi(switches["-do_table_lookup"].c_str())); }

#if TABLE_PRODUCT_QUANTIZE
    if (switches.count("-product_quantize"))  { product_quantize = bool(atoi(switches["-product_quantize"].c_str())); }
    if (switches.count("-product_quantize_log"))  { product_quantize_log = bool(atoi(switches["-product_quantize_log"].c_str())); }
    if (switches.count("-product_quantize_dt_all"))  { product_quantize_dt_all = bool(atoi(switches["-product_quantize_dt_all"].c_str())); }
    if (switches.count("-product_quantize_dims"))  { product_quantize_dims = atoi(switches["-product_quantize_dims"].c_str()); }
    if (switches.count("-product_quantize_mapn"))  { product_quantize_mapn = atoi(switches["-product_quantize_mapn"].c_str()); }
    if (switches.count("-product_quantize_knn"))  { product_quantize_knn = atoi(switches["-product_quantize_knn"].c_str()); }
#endif

#if (TABLE_CLUSTER_KMEANS||TABLE_PRODUCT_QUANTIZE)
    if (switches.count("-kmeans_attempts"))  { kmeans_attempts = atoi(switches["-kmeans_attempts"].c_str()); }
    if (switches.count("-kmeans_max_iters"))  { kmeans_max_iters = atoi(switches["-kmeans_max_iters"].c_str()); }
    if (switches.count("-kmeans_eps"))  { kmeans_eps = atof(switches["-kmeans_eps"].c_str()); }
#endif
    
    if (switches.count("-grid_ndims"))  { grid_ndims = atoi(switches["-grid_ndims"].c_str()); }

    if (switches.count("-table_ratio"))  { table_ratio = atof(switches["-table_ratio"].c_str()); }

    if (switches.count("-prop_dist"))  { prop_dist = atoi(switches["-prop_dist"].c_str()); }

    if (switches.count("-kdtree_add_unique"))  { kdtree_add_unique = bool(atoi(switches["-kdtree_add_unique"].c_str())); }
    
    if (switches.count("-filter_dims")) {
        vector<string> nL = str_split(switches["-filter_dims"], ',');
        ndims = nL.size();
        filter_dims.resize(ndims);
        for (int i = 0; i < (int) nL.size(); i++) {
            filter_dims[i] = atoi(nL[i].c_str());
        }
    }
    
    if (switches.count("-auto_filter")) {
        auto_filter = bool(atoi(switches["-auto_filter"].c_str()));
    }
    
    if (auto_filter) {
        if      (ndims == 6 ) { filter_dims = vector<int>({2,1,5,0,3,10}); }
        else if (ndims == 7 ) { filter_dims = vector<int>({5,1,3,10,0,6,2}); }
        else if (ndims == 8 ) { filter_dims = vector<int>({2,1,6,3,0,5,19,10}); }
        else if (ndims == 10) { filter_dims = vector<int>({2,6,19,10,3,0,1,15,5,13}); }
        else if (ndims == 11) { filter_dims = vector<int>({13,19,6,1,3,30,10,0,5,28,2}); }
        else { fprintf(stderr, "auto_filter -ndims not supported: %d\n", ndims); ASSERT2(false, "ndims not supported"); }
    }
}

PatchTableParams::PatchTableParams() {
    set_defaults();
}

PatchTableParams::PatchTableParams(int argc, const char **argv) {
    set_defaults();
    set_from_switches(argc, argv);
}

PatchTableParams::PatchTableParams(int argc, char **argv) {
    set_defaults();
    set_from_switches(argc, argv);
}

void print_switches() {
    printf("Simplified interface:\n");
    printf("  -speed i                 -- Use preset parameters for speed: 0 (slowest) to 9 (fastest)\n");
    printf("\n");
    printf("Main parameters for controlling speed/accuracy:\n");
    printf("  -ndims d                 -- Number of dimensions for patch descriptor (typically 6-8, can be 20 or more)\n");
    printf("  -grid_ndims d            -- Number of dimensions for grid (typically 6-8)\n");
    printf("  -do_prop b               -- Do propagation (0 or 1)\n");
    printf("  -do_rs b                 -- Do random search (0 or 1)\n");
    printf("\n");
    printf("Other parameters:\n");
    printf("  -nchroma d               -- Use chroma channels (0 means greyscale, 1 means 3-channel)\n");
    printf("  -nslices n1,n2,n3,...    -- Slices for each dimension\n");
    printf("  -ntables n               -- Number of tables\n");
    printf("  -table_knn k             -- Number of k-nearest neighbors for table\n");
    printf("  -run_dt b                -- Whether to run distance transform (0 or 1)\n");
    printf("  -dt_algo raster|prop|brute|kdtree|downsample|hybrid|euclidean -- Distance transform algorithm\n");
    printf("  -dt_iters n              -- Max number of dt iterations (-1 unlimits, default is unlimited)\n");
    printf("  -randomize_dt b          -- Randomize element selected by distance transform (0 or 1)\n");
    printf("  -parallel_dt b           -- Whether to parallelize distance transform (0 or 1)\n");
    printf("  -dt_threads n            -- Thread count for parallel dt (by default uses -threads count)\n");
    printf("  -populate_nearest b      -- Populate nearest grid cell only (0 or 1)\n");
    printf("  -populate_random b       -- Populate in random order (0 or 1)\n");
    printf("  -partition_step n        -- Sample every nth patch in x/y dimensions for building grid cells\n");
    printf("  -coherence_spatial k     -- Values above 0.0 increase spatial coherence of result NNF\n");
    printf("  -coherence_temporal k    -- Values above 0.0 increase temporal coherence of result NNF\n");
    printf("  -recalc_dist_temporal b  -- For temporal coherence, assume patch distances are unreliable; recalculate\n");
    printf("  -colorspace luv|lab|rgb  -- Colorspace to convert to before querying/building table\n");
    printf("  -lookup_algo table|brute|kdtree|treecann -- Lookup algorithm\n");
    printf("  -prop_dist n             -- Additional propagation distance n\n");
    printf("  -check_existing b        -- In match mode checks computed distances (0 or 1)\n");
    printf("  -sanitize_input b        -- Sanitize input for prev_nnf (clamps coords to avoid crashes, 0 or 1)\n");
    printf("  -do_table_lookup b       -- Do lookup in table (0 or 1, default 1)\n");
    printf("  -init_random b           -- Initialize with random NNF (0 or 1, default 0)\n");
    printf("  -prop_iters n            -- Number of query iterations (run propagate/k-coherence algorithm)\n");
    printf("  -triangle_factor f       -- Early termination threshold for triangle inequality (tau**2)\n");
    printf("\n");
#if TABLE_CLUSTER_KMEANS
    printf("  -cluster_kmeans b        -- Use k-means clustering at coarse level and table at fine level\n");
    printf("  -cluster_count n         -- Cluster count for k-means\n");
    printf("  -cluster_limit n         -- Limit number of grid cells within each cluster to n million\n");
    printf("\n");
#endif
#if TABLE_PRODUCT_QUANTIZE
    printf("  -product_quantize b      -- Use product quantization\n");
    printf("  -product_quantize_log b  -- Log product quantization debugging information to files\n");
    printf("  -product_quantize_dims n -- Dimension count for each subspace of product quantization\n");
    printf("  -product_quantize_mapn n -- Slices per dimension for product quantization mapping\n");
    printf("  -product_quantize_knn n  -- knn for product quantization distance transform\n");
    printf("  -product_quantize_dt_all b\n");
    printf("\n");
#endif
#if (TABLE_CLUSTER_KMEANS||TABLE_PRODUCT_QUANTIZE)
    printf("  -kmeans_attempts n\n");
    printf("  -kmeans_max_iters n\n");
    printf("  -kmeans_eps f\n");
    printf("\n");
#endif
    printf("  -treecann_agrid n\n");
    printf("  -treecann_bgrid n\n");
    printf("  -treecann_eps e\n");
    printf("\n");
    printf("  -flann_trees n\n");
    printf("  -flann_checks n\n");
    printf("  -flann_eps e\n");
    printf("  -flann_reset_step n\n");
    printf("  -flann_build_step n\n");
    printf("  -treecann_agrid n\n");
    printf("  -treecann_bgrid n\n");
    printf("\n");
    printf("  -dt_knn k                -- Number of k-NN if dt_algo is kdtree, randomizes neighbor\n");
    printf("  -kcoherence k\n");
    printf("  -kcoherence_min_dist d     -- Min distance from original patch required for k-coherence matches\n");
    printf("  -kcoherence_algo pm|flann|ann -- K-coherence algo\n");
    printf("  -kcoherence_iter n       -- Only run k-coherence on iteration n (-1 to run on all iterations)\n");
    printf("  -kcoherence_step n       -- k-coherence downsampling/step size\n");
    printf("  -kcoherence_enrich b     -- Enrichment in k-coherence (0 or 1)\n");
    printf("  -incoherence b           -- Incoherence for k-coherence (0 or 1)\n");
    printf("  -pm_iters n\n");
    printf("  -pm_enrich b             -- Use enrichment for k-NN PM (0 or 1)\n");
    printf("  -dim_algo wh|pca         -- Dimension reduction algorithm\n");
    printf("  -pca_mean b              -- Use mean for PCA (0 or 1)\n");
    printf("  -kdtree_add_unique b     -- Add only patches in unique bins to kd-tree\n");
#if TABLE_ENABLE_ANN
    printf("\n");
    printf("  -ann b                   -- Use TreeCANN-like algorithm, with ANN library\n");
    printf("  -ann_agrid n             -- ANN a image grid spacing\n");
    printf("  -ann_bgrid n             -- ANN b image grid spacing\n");
    printf("  -ann_eps eps             -- ANN search error eps\n");
#endif
    printf("\n");
    printf("  -query_iters n\n");
    printf("  -pca_samples n\n");
    printf("  -spatial r               -- Spatial search max radius (0 to disable)\n");
    printf("  -query_step n            -- Query is made only every n pixels\n");
//    printf("  -prob_rs p               -- Probability of random search\n");
    printf("  -limit n                 -- Limit number of grid cells to n million\n");
    printf("  -filter_dims d1,d2,...   -- Use selected dimensions (also sets ndims to length of the list)\n");
    printf("  -auto_filter b           -- Automatically set filter_dims (0 or 1)\n");
    printf("\n");
    printf("  -threads n               -- OpenMP thread count\n");
    printf("  -load filename           -- Load table from filename\n");
    printf("  -save filename           -- Save table to filename\n");
    printf("  -compare filename.pfm    -- Compare matches with NNF file\n");
}

/* ------------------------------------------------------------------------------------------------
   Serialization of basic data-types
   ------------------------------------------------------------------------------------------------ */

#define LOAD_SAVE_DTYPE(dtype) \
void save(FILE *f, dtype v) { \
    ASSERT2(fwrite((void *) &v, sizeof(dtype), 1, f) == 1, "expected to write 1 element"); \
} \
\
void load(FILE *f, dtype &v) { \
    ASSERT2(fread((void *) &v, sizeof(dtype), 1, f) == 1, "expected to read 1 element"); \
} \
\
void save(FILE *f, const vector<dtype> &v) { \
    if (TABLE_DEBUG_SERIALIZE) { printf("save(vector): saving %d elements each of size %d\n", v.size(), sizeof(dtype)); } \
    save(f, int(v.size())); \
    ASSERT2(fwrite((void *) &v[0], sizeof(dtype), v.size(), f) == v.size(), "save(vector): expected to write n elements"); \
} \
\
void load(FILE *f, vector<dtype> &v) { \
    int size; \
    load(f, size); \
    v.resize(size); \
    if (TABLE_DEBUG_SERIALIZE) { printf("load(vector): loading %d elements each of size %d\n", v.size(), sizeof(dtype)); } \
    ASSERT2(fread((void *) &v[0], sizeof(dtype), size, f) == size, "load(vector): expected to read n elements"); \
} \
\
void save(FILE *f, const Array<dtype> &v) { \
    save(f, v.sizes); \
    if (TABLE_DEBUG_SERIALIZE) { printf("save(Array): saving %d elements each of size %d\n", v.nelems, sizeof(dtype)); } \
    ASSERT2(fwrite((void *) &v.data[0], sizeof(dtype), v.nelems, f) == v.nelems, "save(Array): expected to write nelems elements"); \
} \
\
void load(FILE *f, Array<dtype> &v) { \
    vector<int> sizes; \
    load(f, sizes); \
    v.resize(sizes); \
    if (TABLE_DEBUG_SERIALIZE) { printf("load(Array): loading %d elements each of size %d\n", v.nelems, sizeof(dtype)); } \
    ASSERT2(fread((void *) &v.data[0], sizeof(dtype), v.nelems, f) == v.nelems, "load(Array): expected to read nelems elements"); \
}

LOAD_SAVE_DTYPE(int32_t);
LOAD_SAVE_DTYPE(char);
LOAD_SAVE_DTYPE(int8_t);
LOAD_SAVE_DTYPE(int16_t);
LOAD_SAVE_DTYPE(int64_t);
LOAD_SAVE_DTYPE(float);
LOAD_SAVE_DTYPE(double);

void save(FILE *f, const string &s) {
    save(f, int(s.size()));
    for (int i = 0; i < (int) s.size(); i++) {
        save(f, char(s[i]));
    }
}

void save(FILE *f, const char *s) {
    save(f, string(s));
}

void load(FILE *f, string &s) {
    int size;
    load(f, size);
    if (TABLE_DEBUG_SERIALIZE) { printf("load(string): length is %d\n", size); }
    char *buf = new char[size];
    for (int i = 0; i < size; i++) {
        load(f, buf[i]);
    }
    s = string(buf, size);
}

void load_id(FILE *f, const string &s) {
    string s2;
    load(f, s2);
    if (s2 != s) {
        fprintf(stderr, "Error: load_id expected |%s|, received |%s|\n", s.c_str(), s2.c_str());
        ASSERT2(false, "load from file failed: load_id expected matching string");
    }
}

AdjacencySet::AdjacencySet(int n, bool unique_) :sets(unique_ ? 0: n, NULL), setsL(n, NULL), unique(unique_) {
}

AdjacencySet::~AdjacencySet() {
    for (int i = 0; i < (int) sets.size(); i++) {
        delete sets[i];
    }
    for (int i = 0; i < (int) setsL.size(); i++) {
        delete setsL[i];
    }
}
    
void AdjacencySet::add(int i, int j) {
    if (!unique) {
        if (sets[i] == NULL) {
            sets[i] = new unordered_set<int>();
        }
        sets[i]->insert(j);
    } else {
        if (setsL[i] == NULL) {
            setsL[i] = new vector<int>();
        }
        setsL[i]->push_back(j);
    }
}

void AdjacencySet::compute_sets() {
    if (!unique) {
        for (int i = 0; i < (int) sets.size(); i++) {
            if (sets[i] && !setsL[i]) {
                setsL[i] = new vector<int>(sets[i]->begin(), sets[i]->end());
            }
        }
    }
}

vector<int> scale_slices_to_limit(vector<int> nslices, double limit) {
    if (limit > 0) {
        vector<int> nslices0(nslices);
        double factor = 1.0;
        for (int i = 0; i < (int) nslices.size(); i++) {
            nslices[i] = MAX(nslices[i], 1);
        }
        while (product(nslices) > limit) {
            factor *= 0.995;
            nslices = multiply(nslices0, factor);
            for (int i = 0; i < (int) nslices.size(); i++) {
                nslices[i] = MAX(nslices[i], 1);
            }
        }
        for (;;) {
            nslices0 = nslices;
            bool success = false;
            for (int i = 0; i < (int) nslices.size(); i++) {
                nslices[i]++;
                if (product(nslices) > limit) {
                    nslices = nslices0;
                } else {
                    success = true;
                    nslices0 = nslices;
                }
            }
            if (!success) { break; }
        }
//        printf("got -limit switch, factor=%f, product(nslices) = %f\n", factor, product(nslices));
    }
    return nslices;
}
