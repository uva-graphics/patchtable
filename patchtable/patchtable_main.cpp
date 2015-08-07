
/* ------------------------------------------------------------------
   Main program to test PatchTable and GCK
   ------------------------------------------------------------------

   - Video analogies
     - Higher priority than optimizations

   - Optimizations:
     - PM:
        - Vanilla mode PM (no special if branches)
        - Compute correlation instead of distance
        - Pad WH vectors to multiple of 4 so masking is not needed in patch_dist_approx()
     - Downsample/upsample algorithm for distance transform
     - Profiling
     - Benchmark query 1 vs N images. Precomputation time, query time, total time.

   - Applications:
     - Inpainting with billions of patches
       - Maybe "repainting" not inpainting -- redesigning or inpainting a selected region
       - Interactivity: lighter/darker and color sliders, change season (draw snow), draw color in the sky, draw autumn leaves, draw edges
     - Inpainting / edit propagation in personal photo collections -- NRDC, propagate inpaintings
     - NRDC 
     - Video?
        - Fast optical flow? But optical flow has a lot of competition. So maybe fast NRDC or matching between stereo or light field video frames.
           - For search windows could add x/y as dimensions to the table. But for small search windows brute force is probably good enough.
        - NPR with video analogies?
     - Fast synthesis with rotations+scales?
     
   - Optimizations:
     - Use nth_element for grid spacing instead of sort? Try also increasing -partition_step on the benchmark
     - Parallel dt, and parallel GCK
     - Find a way to avoid the problem where a given bin is never used (due to many patch WH values repeating).
     - Rotations+scales?
     - Follow on paper: Can consider irregular grid to be a nonlinear distortion applied to each dimension independently.
       Could also apply to m=2 or more dimensions simultaneously a nonlinear distortion, via an m-dimensional lookup table.
*/
#include "patchtable.h"
#include "analogy.h"

/* ------------------------------------------------------------------------------------------------
   PatchTable test functions
   ------------------------------------------------------------------------------------------------ */

#define TABLE_TEST_NORMAL_GCK 1

template<class in_type, class out_type, int n>
void test_gck(int argc, char **argv, string msg, bool save=true) {
    Array<in_type> a;
    a = load_color_image<in_type>(argv[1]);
    Array<out_type> b;
    for (int i = 0; i < 5; i++) {
        double T0 = wall_time();
#if TABLE_TEST_NORMAL_GCK
        gck_8<in_type, out_type, n, n>(a, b);
#else
        IncrementalGCK_8<in_type, out_type, n, n>(a, b).compute_all();
#endif
        double T1 = wall_time();
        printf("gck (%s, %d): %f secs\n", msg.c_str(), n, T1-T0);
    }
    if (save) {
        save_color_image(b, argv[2]);
    }
    printf("\n");
}

void usage() {
    printf("patchtable gck in.png out.pfm\n");
    printf("patchtable match a.png b.png out.pfm [-prev_nnf prev.pfm] [-allowed_patches mask.png] [-is_descriptor 0|1]\n");
    printf("patchtable analogy a.png aprime.png b.png|b%%05d.png out.png|out%%05d.png\n");
    printf("  If %%d format is used, processes video starting at frame 1\n");
    printf("patchtable test_descriptor\n");
    printf("patchtable check_dist a.png b.png nnf.pfm [nnf_dist_recalc.pfm] -- Prints average distance for a -> b nnf\n");
    printf("\n");
    print_switches();
    print_analogy_switches();
    printf("\n");
    printf("Examples:\n");
    printf("\n");
    printf("  Match with accuracy roughly similar to PatchMatch:\n");
    printf("  $ ./patchtable match vidpair0/a.png vidpair0/b.png out.pfm\n");
    printf("\n");
    printf("  Match with increased coherence (mimics smooth NNFs of PatchMatch):\n");
    printf("  $ ./patchtable match vidpair0/a.png vidpair0/b.png out.pfm -speed 3 -coherence_spatial 4\n");
    printf("\n");
    exit(1);
}

double get_mean_dist(const Array<double> &ann) {
    double mean_dist = 0;
    for (int y = 0; y < ann.height(); y++) {
        for (int x = 0; x < ann.width(); x++) {
            double d = double(ann(y, x, NNF_DIST));
            if (d < 0) { d = 0; }
            mean_dist += sqrt(d);
        }
    }
    mean_dist = mean_dist/double(ann.height()*ann.width());
    return mean_dist;
}

template<int is_descriptor, int check_existing=0>
double recompute_mean_dist(const Array<float> &a, const Array<float> &b, const Array<double> &ann, PatchTableParams *p, bool overwrite=false) {
    int ann_w = a.width()-p->patch_w+1;
    int ann_h = a.height()-p->patch_w+1;
    int bew = b.width()-p->patch_w+1;
    int beh = b.height()-p->patch_w+1;
    if (is_descriptor) {
        ann_w = a.width();
        ann_h = a.height();
        bew = b.width();
        beh = b.height();
    }
    if (ann.height() != ann_h || ann.width() != ann_w) {
        fprintf(stderr, "recompute_mean_dist: ann_w=%d, ann_h=%d, ann.width()=%d, ann.height()=%d, a.width()=%d, a.height()=%d, is_descriptor=%d\n", ann_w, ann_h, ann.width(), ann.height(), a.width(), a.height(), is_descriptor);
        ASSERT2(false, "expected ann size to match a size (after removing patch_w-1 from a size)");
    }
    
    double sum_dist = 0.0;
    
    for (int ay = 0; ay < ann_h; ay++) {
        for (int ax = 0; ax < ann_w; ax++) {
            int bx = ann(ay, ax, NNF_X);
            int by = ann(ay, ax, NNF_Y);
//            if (!in_bounds(bx, bew) || !in_bounds(by, beh)) { fprintf(stderr, "nnf out of bounds: %d, %d => %d, %d, b size is %dx%d\n", ax, ay, bx, by, b.width(), b.height()); exit(1); }
            if (bx < 0) { bx = 0; }
            else if (bx >= bew) { bx = bew-1; }
            if (by < 0) { by = 0; }
            else if (by >= beh) { by = beh-1; }
            double dist = 0.0;
            if (!is_descriptor) {
                for (int dy = 0; dy < p->patch_w; dy++) {
                    for (int dx = 0; dx < p->patch_w; dx++) {
                        for (int channel = 0; channel < a.channels(); channel++) {
                            ASSERT(in_bounds(ay+dy, a.height()), "ay+dy out of bounds a.height()");
                            ASSERT(in_bounds(ax+dx, a.width()),  "ax+dx out of bounds a.width()");
                            ASSERT(in_bounds(by+dy, b.height()), "by+dy out of bounds b.height()");
                            ASSERT(in_bounds(bx+dx, b.width()),  "bx+dx out of bounds b.width()");
                            ASSERT(in_bounds(channel, b.channels()), "channel out of bounds b.channels()");
                            double delta = a(ay+dy, ax+dx, channel) - b(by+dy, bx+dx, channel);
    //                        if (ax % 50 == 0 && ay == 0) {
    //                            printf("%d %d: %f %f, delta=%f\n", dx, dy, a(ay+dy, ax+dx, channel), b(by+dy, bx+dx, channel), delta);
    //                        }
                            dist += delta*delta;
                        }
                    }
                }
            } else {
                for (int i = 0; i < a.channels(); i++) {
                    double delta = a(ay, ax, i) - b(by, bx, i);
                    dist += delta*delta;
                }
            }
//            if (ax % 50 == 0 && ay == 0) {
//                printf("d(%d, %d) = %f, bx, by=%d, %d\n\n", ax, ay, sqrt(dist), bx, by);
//            }
//            ann(ay, ax, NNF_DIST) = dist;
            if (overwrite) {
                ann.get_nearest(ay, ax, NNF_DIST) = dist;
            }
            if (check_existing || p->check_existing) {
                double d_existing = ann(ay, ax, NNF_DIST);
                if (fabs(dist - d_existing) > 1e-4) { fprintf(stderr, "recomputed distance %f differs from existing distance %f at %d, %d\n", dist, d_existing, ax, ay); ASSERT2(false, "error"); }
            }
            
            sum_dist += sqrt(dist);
        }
    }
    
    return sum_dist/double(ann_w*ann_h);
}

void test_descriptor(int w=400, int h=300, int k=8) {
    Array<float> a0, b0;
    vector<int> size({h, w, k});
    a0.assign(Array<float>::random(size));
    b0.assign(Array<float>::random(size));
    
    Array<float> a, b;
    a.assign(a0);
    b.assign(b0);
    
    PatchTableParams *p = new PatchTableParams();
    p->is_descriptor = true;
    vector<double> dL(3);
    
    for (int reduce = 0; reduce < 3; reduce++) {
        if (reduce == 1) {
            p->grid_ndims = k-2;
        }
        else if (reduce == 2) {
            p->grid_ndims = p->ndims = k-2;
            Array<float> ap, bp;
            vector<int> start({0, 0, 0});
            vector<int> extent({h, w, k-2});
            ap.resize(extent);
            bp.resize(extent);
            a.copy_rect(ap, start, start, extent);
            b.copy_rect(bp, start, start, extent);
            a.assign(ap);
            b.assign(bp);
        }
        printf("test_descriptor: reduce=%d, ndims=%d, b channels: %d\n", reduce, p->grid_ndims, b.channels());
        
        PatchTable<> btable(p, b);
        printf("test_descriptor: done constructing table\n"); fflush(stdout);
        Array<double> ann;
        btable.lookup(a, ann);
        
        if (ann.width() != a.width() || ann.height() != a.height()) {
            fprintf(stderr, "expected ann size %dx%d to match a size %dx%d\n", ann.width(), ann.height(), a.width(), a.height()); exit(1);
        }
        for (int y = 0; y < ann.height(); y++) {
            for (int x = 0; x < ann.width(); x++) {
                int bx(ann(y, x, NNF_X));
                int by(ann(y, x, NNF_Y));
                if ((unsigned) bx >= (unsigned) b.width() || (unsigned) by >= (unsigned) b.height()) {
                    fprintf(stderr, "expected ax, ay=%d, %d => bx, by=%d, %d to be in bounds %dx%d\n", x, y, bx, by, b.width(), b.height());
                    exit(1);
                }
            }
        }

        dL[reduce] = recompute_mean_dist<1>(a0, b0, ann, p); //get_mean_dist(ann);
        printf("mean_dist: %f (reduce=%d)\n", dL[reduce], reduce);
    }
    if (dL[0] < dL[1] && dL[1] < dL[2] && dL[0] < 0.77) {
        printf("test_descriptor: OK\n");
    } else {
        fprintf(stderr, "test_descriptor: Failed, distances not as expected\n"); exit(1);
    }
}

void test_match(const char *a_filename, const char *b_filename, const char *out_filename, int argc, char **argv, double image_scale) {
    PatchTableParams *p = new PatchTableParams(argc, argv);
    p->calc_exact_dist = true;
    
    map<string, string> switches = parse_switches(argc, argv);
    bool is_descriptor = false;
    if (switches.count("-is_descriptor")) { is_descriptor = bool(atoi(switches["-is_descriptor"].c_str())); }
    p->is_descriptor = is_descriptor;
    
    Array<double> *prev_nnf = NULL;
    if (switches.count("-prev_nnf")) {
        prev_nnf = new Array<double>(load_color_image<double>(switches["-prev_nnf"].c_str()));
    }
    Array<int32_t> *allowed_patches = NULL;
    if (switches.count("-allowed_patches")) {
        Array<double> allowed_patches_f(load_color_image<double>(switches["-allowed_patches"].c_str()));
        ASSERT2(allowed_patches_f.dimensions() >= 3, "expected allowed_patches image to have >= 3 dimensions");
        allowed_patches = new Array<int32_t>(allowed_patches_f.height(), allowed_patches_f.width());
        for (int y = 0; y < allowed_patches_f.height(); y++) {
            for (int x = 0; x < allowed_patches_f.width(); x++) {
                (*allowed_patches)(y, x) = allowed_patches_f(y, x, 0);
            }
        }
    }
    
    if (p->verbose) {
        printf("----------------------------------------\n");
        printf("Precomputation:\n");
        printf("----------------------------------------\n");
        printf("\n");
    }
    Array<float> a(load_color_image<float>(a_filename));
    Array<float> b(load_color_image<float>(b_filename));
    if (p->verbose) {
        printf("a size: width %d x height %d x channels %d, %s\n", a.width(), a.height(), a.channels(), a_filename);
        printf("b size: width %d x height %d x channels %d, %s\n", b.width(), b.height(), b.channels(), b_filename);
    }
    Array<float> a_desc, b_desc;
    if (is_descriptor) {
        if (a.channels() == 3) {
            gck_no_pad<float, float>(a, a_desc, p->ndims-p->nchroma*2, p->nchroma, p->patch_w);
            gck_no_pad<float, float>(b, b_desc, p->ndims-p->nchroma*2, p->nchroma, p->patch_w);
        } else {
            a_desc.assign(a);
            b_desc.assign(b);
        }
    }

    double T0_table = wall_time();
    PatchTable<> table(p, is_descriptor ? b_desc: b, allowed_patches);
    double T1_table = wall_time();
//    if (p->verbose && p->is_descriptor) {
//        printf("a_desc size: %s, b_desc size: %s\n", vector_to_str_int(a_desc.sizes).c_str(), vector_to_str_int(b_desc.sizes).c_str());
//    }
    Array<double> ann;
    double latest_time = 1e100;
    for (int iter = 0; iter < p->query_iters; iter++) {
        if (p->verbose) {
            printf("----------------------------------------\n");
            printf("Query #%d\n", iter);
            printf("----------------------------------------\n");
            printf("\n");
        }
        
        latest_time = table.lookup(is_descriptor ? a_desc: a, ann, prev_nnf);
        
        /* Explanation of ann (nearest neighbor field matching from a -> b):
           int b_patch_x       = ann(a_patch_y, a_patch_x, NNF_X);
           int b_patch_y       = ann(a_patch_y, a_patch_x, NNF_Y);
           double b_patch_dist = ann(a_patch_y, a_patch_x, NNF_DIST); */

        if (p->verbose) {
            printf("\n");
        }
    }

    save_color_image(ann, out_filename);

    double total_time = latest_time+(T1_table-T0_table);
    double mean_dist = get_mean_dist(ann);
    double mean_dist_recomputed = is_descriptor ? recompute_mean_dist<1>(a_desc, b_desc, ann, p): recompute_mean_dist<0>(a, b, ann, p);
    if (p->is_descriptor) { mean_dist = mean_dist_recomputed; }
    if (p->verbose) {
        printf("Combined precompute and lookup time: %f\n", total_time);
        printf("mean_dist: %f\n", mean_dist);
        printf("mean_dist_recomputed: %f\n", mean_dist_recomputed);
    } else {
        printf("total_time: %f, query_time: %f, mean_dist: %f, mean_dist_recomputed: %f\n", total_time, latest_time, mean_dist, mean_dist_recomputed);
    }

	delete p;
    if (prev_nnf) { delete prev_nnf; }
    if (allowed_patches) { delete allowed_patches; }
}

void check_dist(const char *a_filename, const char *b_filename, const char *nnf_filename, int argc, char **argv) {
    Array<float> a(load_color_image<float>(a_filename));
    Array<float> b(load_color_image<float>(b_filename));
    Array<double> ann(load_color_image<double>(nnf_filename));
    ASSERT2(ann.channels() == 3, "expected ann to have 3 channels");

    PatchTableParams *p = new PatchTableParams(argc, argv);

    bool check = false;
    map<string, string> sw = parse_switches(argc, argv);
    if (sw.count("-check")) {
        check = bool(atoi(sw["-check"].c_str()));
    }
    bool overwrite = argc > 4;
    if (check) {
        printf("%f\n", recompute_mean_dist<0, 1>(a, b, ann, p, overwrite));
    } else {
        printf("%f\n", recompute_mean_dist<0, 0>(a, b, ann, p, overwrite));
    }
    if (overwrite) {
        char *out_filename = argv[4];
        save_color_image<double>(ann, out_filename);
    }
}

/* ------------------------------------------------------------------------------------------------
   Main program
   ------------------------------------------------------------------------------------------------ */

int main(int argc0, char **argv0) {
    int argc = argc0-1;
    char **argv = argv0+1;

    if (argc < 1) {
        usage();
    }
    
    if (!strcmp(argv[0], "gck")) {
        ASSERT2(argc >= 3, "expected 3 args in gck mode");
        for (int i = 0; i < 2; i++) {
            if (i == 0) {
                test_gck<double,  double,  10>(argc, argv,  "double, double",   false);
                test_gck<float,   float,   10>(argc, argv,  "float, float",     false);
                test_gck<uint8_t, int32_t, 10>(argc, argv,  "uint8_t, int32_t", false);
            } else {
                test_gck<double,  double,  64>(argc, argv,  "double, double",   false);
                test_gck<float,   float,   64>(argc, argv,  "float, float",     false);
                test_gck<uint8_t, int32_t, 64>(argc, argv,  "uint8_t, int32_t", true);
            }
        }
    } else if (!strcmp(argv[0], "match")) {
        ASSERT2(argc >= 4, "expected 4 args in match mode");
        test_match(argv[1], argv[2], argv[3], argc, argv, 1.0/255);
    } else if (!strcmp(argv[0], "analogy")) {
        analogy(argc-1, argv+1);
    } else if (!strcmp(argv[0], "test_descriptor")) {
        test_descriptor();
    } else if (!strcmp(argv[0], "check_dist")) {
        ASSERT2(argc > 3, "expected 3 args in check_dist mode");
        check_dist(argv[1], argv[2], argv[3], argc, argv);
    } else {
        usage();
    }
    
    return 0;
}
