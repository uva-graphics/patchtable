
#ifndef _patchtable_params_h
#define _patchtable_params_h

#include <vector>
#include <string>

using std::vector;
using std::string;

#define TABLE_PRODUCT_QUANTIZE  0           /* Enable product quantization */
#define TABLE_RS_TABLE          0
#define TABLE_CLUSTER_KMEANS    0           /* Use k-means clustering as a coarse level and tables inside fine level */

#define gck_ymin(a) (p->patch_w-1)
#define gck_xmin(a) (p->patch_w-1)
#define gck_xmax(a) ((a).width()-p->patch_w+1)
#define gck_ymax(a) ((a).height()-p->patch_w+1)
#define gck_ew(a) (gck_xmax(a) - gck_xmin(a))
#define gck_eh(a) (gck_ymax(a) - gck_ymin(a))

/* Run patchtable program, or see patchtable.cpp function print_switches() for documentation of the parameters */
class PatchTableParams { public:
    bool randomize_dt;
    bool parallel_dt;
    double limit;
    int threads;
    int dt_threads;
    bool convert_colorspace;
    int colorspace;                           // One of PATCHTABLE_COLORSPACE_*.
    bool regular_grid;                        // No longer used (deprecated)
    int verbose;
    int table_knn;
    bool populate_nearest;
    bool populate_random;
    int partition_step;
    bool check_existing;
    bool do_table_lookup;
    bool init_random;
    double treecann_eps;
    int allowed_index;                        // When allowed_patches(y, x) == allowed_index, a patch is allowed
    
#if TABLE_CLUSTER_KMEANS
    bool cluster_kmeans;
    int cluster_count;
    double cluster_limit;
#endif

#if TABLE_PRODUCT_QUANTIZE
    bool product_quantize;
    bool product_quantize_log;
    int product_quantize_dims;
    int product_quantize_mapn;
    int product_quantize_knn;
    bool product_quantize_dt_all;
#endif

#if (TABLE_CLUSTER_KMEANS||TABLE_PRODUCT_QUANTIZE)
    int kmeans_attempts;
    int kmeans_max_iters;
    double kmeans_eps;
#endif
    
    int ndims;
    int nchroma;
    int ntables;
    vector<int> nslices;
    bool populate_exptime;
    int dt_mode;
    int dt_algo;
    int dt_iters;
    bool run_dt;
    int query_step;
    int grid_ndims;                           // If this is not -1 (the default) and is_descriptor is true, then uses grid_ndims < ndims
                                              // to build the table, but compares distances using the original ndims.
    
    int prop_dist;
    bool calc_exact_dist;
    bool do_prop;
    int prop_y_step;
    bool do_rs;
    int spatial;
    double prob_rs;                           // No longer used (deprecated)
    bool do_calc_dist;
    int lookup_algo;
    
    bool is_descriptor;                       // If true assumes input is a h x w x ndims feature descriptor that is matched directly
                                              // (no extraction of patches or dimension reduction is performed)
    bool descriptor_padded;                   // If true assumes input is padded on all sides the same as gck.h
    bool sanitize_input;                      // If true then cleans up prev_nnf input by clamping it if necessary
    
    string load_filename;
    string save_filename;
    
    int patch_w;
    double coherence_spatial;                 // Spatial coherence -- Values above 0.0 increase coherency of output NNF
    double coherence_temporal;                // Temporal coherence -- Values above 0.0 increase temporal coherence, if ann_prev provided to lookup()
    bool recalc_dist_temporal;                // When using temporal coherence assume input patch distances are unreliable and recalculate them
    
    double table_ratio;
    
    int dt_knn;
    
    int treecann_agrid;
    int treecann_bgrid;
    
    int flann_build_step;
    int flann_trees;
    int flann_checks;
    double flann_eps;
    int flann_reset_step;
    
    int query_iters;
    
    int prop_iters;
    double enrich_limit;
    
    int kcoherence;
    int kcoherence_min_dist;
    int kcoherence_algo;
    int kcoherence_iter;
    int kcoherence_step;
    int kcoherence_enrich;
    bool save_kcoherence;
    bool incoherence;
    int incoherence_min_dist;
    double triangle_factor;
    
    int pm_iters;                             // Separate parameters for core PM
    int pm_knn;
    int pm_min_dist;
    bool pm_enrich;
    
    int pca_samples;
    bool pca_mean;
    bool pca_subtract_mean;
    int dim_algo;
    vector<int> filter_dims;
    bool auto_filter;
    string compare;
    
    bool kdtree_add_unique;

    bool ann;
    int ann_agrid;
    int ann_bgrid;
    double ann_eps;

#if TABLE_RS_TABLE
    int rs_table;
#endif

    void set_speed(int i);
    void set_defaults();
    void set_from_switches(int argc, char **argv);
    void set_from_switches(int argc, const char **argv);
    void set_from_switches(string s);                   // Parses switches from string, e.g. "-limit 1 -do_rs 1"
    PatchTableParams();                       /* Sets parameters to defaults */
    PatchTableParams(int argc, char **argv);        /* Sets parameters based on command-line switches */
    PatchTableParams(int argc, const char **argv);  /* Sets parameters based on command-line switches */
};

#endif

