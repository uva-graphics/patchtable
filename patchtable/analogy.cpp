
/* TODO:
    - Check if lowering speed improves synthesis quality.
    - Optimization: Could only process 1 channel when building pyramids / upsampling.
    - Optimization: Use only 1 correction iteration at finest level(s), and use faster parameters at finest level(s).
    - Tune min_size and optical flow parameters to get best quality/fastest results. */

#include "analogy.h"
#include "patchtable.h"


#define ANALOGY_DEFAULT_PATCH_W      8
#define ANALOGY_MAX_FLOW_DIST        0          /* Support max_flow_dist parameter */

#define ANALOGY_ERROR_OUT_OF_BOUNDS  0

#define OPTICAL_FLOW_SIMPLEFLOW      0
#define OPTICAL_FLOW_DUALTV          1
#define OPTICAL_FLOW_FARNEBACK       2

class AnalogyParams { public:
    int lock;
    int min_size;
    int levels;
    int correction_iters;
    bool log;
    bool vote_overlapping;
    bool lum_only;
    bool analogy_verbose;
    double ap_weight;
    
#if ANALOGY_MAX_FLOW_DIST
    double max_flow_dist;
#endif
    
    bool optical_flow;
    int of_algo;

    double of_tau;
    double of_lambda;
    double of_theta;
    int of_nscales;
    int of_warps;
    double of_epsilon;
    int of_iterations;

    int of_layers;
    int of_block_size;
    int of_max_flow;

    double of_pyr_scale;
    int of_levels;
    int of_winsize;
    int of_poly_n;
    double of_poly_sigma;

    double of_downsample;
    double upweight_coherence;

    int update_every;
    int start_index;
    
    int fine_levels;
    double fine_coherence_spatial;
    double fine_upweight_coherence;

    double shot_detect;
    
    void set_from_switches(int argc, char **argv) {
        map<string, string> switches = parse_switches(argc, argv);

        if (switches.count("-lock"))                { lock = atoi(switches["-lock"].c_str()); }
        if (switches.count("-shot_detect"))         { shot_detect = atof(switches["-shot_detect"].c_str()); }
        if (switches.count("-min_size"))            { min_size = atoi(switches["-min_size"].c_str()); }
        if (switches.count("-levels"))              { levels = atoi(switches["-levels"].c_str()); }
        if (switches.count("-correction_iters"))    { correction_iters = atoi(switches["-correction_iters"].c_str()); }
        if (switches.count("-log"))                 { log = bool(atoi(switches["-log"].c_str())); }
        if (switches.count("-vote_overlapping"))    { vote_overlapping = bool(atoi(switches["-vote_overlapping"].c_str())); }
        if (switches.count("-lum_only"))            { lum_only = bool(atoi(switches["-lum_only"].c_str())); }
        if (switches.count("-analogy_verbose"))     { analogy_verbose = bool(atoi(switches["-analogy_verbose"].c_str())); }
        if (switches.count("-upweight_coherence"))  { upweight_coherence = atof(switches["-upweight_coherence"].c_str()); }
        if (switches.count("-ap_weight"))           { ap_weight = atof(switches["-ap_weight"].c_str()); }

        if (switches.count("-fine_levels"))             { fine_levels = atoi(switches["-fine_levels"].c_str()); }
        if (switches.count("-fine_coherence_spatial"))  { fine_coherence_spatial = atof(switches["-fine_coherence_spatial"].c_str()); }
        if (switches.count("-fine_upweight_coherence")) { fine_upweight_coherence = atof(switches["-fine_upweight_coherence"].c_str()); }

#if ANALOGY_MAX_FLOW_DIST
        if (switches.count("-max_flow_dist"))       { max_flow_dist = atof(switches["-max_flow_dist"].c_str()); }
#endif
        
        if (switches.count("-update_every"))     { update_every = atoi(switches["-update_every"].c_str()); }
        if (switches.count("-start_index"))      { start_index = atoi(switches["-start_index"].c_str()); }

        if (switches.count("-optical_flow"))     { optical_flow = bool(atoi(switches["-optical_flow"].c_str())); }
        
        if (switches.count("-of_tau"))           { of_tau = atof(switches["-of_tau"].c_str()); }
        if (switches.count("-of_lambda"))        { of_lambda = atof(switches["-of_lambda"].c_str()); }
        if (switches.count("-of_theta"))         { of_theta = atof(switches["-of_theta"].c_str()); }
        if (switches.count("-of_nscales"))       { of_nscales = atoi(switches["-of_nscales"].c_str()); }
        if (switches.count("-of_warps"))         { of_warps = atoi(switches["-of_warps"].c_str()); }
        if (switches.count("-of_epsilon"))       { of_epsilon = atof(switches["-of_epsilon"].c_str()); }
        if (switches.count("-of_iterations"))    { of_iterations = atoi(switches["-of_iterations"].c_str()); }

        if (switches.count("-of_layers"))        { of_layers = atoi(switches["-of_layers"].c_str()); }
        if (switches.count("-of_block_size"))    { of_block_size = atoi(switches["-of_block_size"].c_str()); }
        if (switches.count("-of_max_flow"))      { of_max_flow = atoi(switches["-of_max_flow"].c_str()); }

        if (switches.count("-of_pyr_scale"))     { of_pyr_scale = atof(switches["-of_pyr_scale"].c_str()); }
        if (switches.count("-of_levels"))        { of_levels = atoi(switches["-of_levels"].c_str()); }
        if (switches.count("-of_winsize"))       { of_winsize = atoi(switches["-of_winsize"].c_str()); }
        if (switches.count("-of_poly_n"))        { of_poly_n = atoi(switches["-of_poly_n"].c_str()); }
        if (switches.count("-of_poly_sigma"))    { of_poly_sigma = atof(switches["-of_poly_sigma"].c_str()); }

        if (switches.count("-of_downsample"))    { of_downsample = atof(switches["-of_downsample"].c_str()); }
        
        if (switches.count("-of_algo")) {
            if (switches["-of_algo"] == string("simpleflow")) { of_algo = OPTICAL_FLOW_SIMPLEFLOW; }
            else if (switches["-of_algo"] == string("dualtv")) { of_algo = OPTICAL_FLOW_DUALTV; }
            else if (switches["-of_algo"] == string("farneback")) { of_algo = OPTICAL_FLOW_FARNEBACK; }
            else { fprintf(stderr, "Unknown -of_algo %s\n", switches["-of_algo"].c_str()); exit(1); }
        }
    }
    
    void set_defaults() {
        lock = 0;
        shot_detect = 0.05;
        
        min_size = 15;
        levels = 8;
        correction_iters = 2;
        log = false;
        vote_overlapping = true;
        lum_only = true;
        analogy_verbose = false;
        upweight_coherence = 0.0;
        start_index = 1;
        ap_weight = 1.0;
#if ANALOGY_MAX_FLOW_DIST
        max_flow_dist = 50.0;
#endif
        
        fine_levels = 2;
        fine_coherence_spatial = -1;
        fine_upweight_coherence = -1;

        optical_flow = false;
        of_algo = OPTICAL_FLOW_FARNEBACK;
        
        update_every = 1;
        
        of_tau = 0.25;
        of_lambda = 0.15;
        of_theta = 0.3;
        of_nscales = 5;
        of_warps = 5;
        of_epsilon = 0.01;
        of_iterations = -1;
        
        of_layers = 3;
        of_block_size = 2;
        of_max_flow = 4;
        
        of_pyr_scale = 0.5;
        of_levels = 3;
        of_winsize = 15;
        of_poly_n = 5;
        of_poly_sigma = 1.2;
        
        of_downsample = 1.0;
    }
    
    AnalogyParams() {
        set_defaults();
    }
    
    AnalogyParams(int argc, char **argv) {
        set_defaults();
        set_from_switches(argc, argv);
    }
};

void print_analogy_switches() {
    printf("\n");
    printf("Analogy parameters:\n");
    printf("  -lock n                    -- Lock synthesis for coarsest n scales if there is optical flow correspondence\n");
    printf("  -shot_detect f             -- If f > 0, sets shot detection threshold\n");
    printf("  -min_size n                -- Minimum size in pyramid\n");
    printf("  -levels n                  -- Pyramid levels\n");
    printf("  -correction_iters n        -- Correction iterations\n");
    printf("  -log b                     -- Log images to disk (0 or 1)\n");
    printf("  -vote_overlapping b        -- Vote using overlapping patches (0 or 1)\n");
    printf("  -lum_only b                -- Synthesize luminance only (0 or 1)\n");
    printf("  -analogy_verbose b         -- Whether analogy should be verbose (0 or 1)\n");
    printf("  -update_every n            -- Update NNF every n frames\n");
    printf("  -ap_weight w               -- Weight on ap (a prime) image in descriptor\n");
    printf("  -upweight_coherence p      -- Upweight coherence by some power p\n");
    printf("\n");
    printf("  -fine_levels n             -- Number of levels considered to be 'fine' scales\n");
    printf("  -fine_coherence_spatial k  -- Spatial coherence at fine scales\n");
    printf("  -fine_upweight_coherence p -- Spatial coherence at fine scales\n");
#if ANALOGY_MAX_FLOW_DIST
    printf("  -max_flow_dist d           -- Decrease temporal coherence to 0 after d pixels of optical flow\n");
#endif
    printf("\n");
    printf("  -optical_flow b            -- Correct NNF with optical flow (0 or 1)\n");
    printf("  -of_algo simpleflow|dualtv|farneback -- OF algorithm\n");
    printf("  -of_tau f, -of_lambda f, -of_theta f, -of_nscales n, -of_warps n, -of_epsilon f, -of_iterations n\n");
    printf("  -of_layers n, -of_block_size n, -of_max_flow n\n");
    printf("  -of_pyr_scale f, -of_levels n, -of_winsize n, -of_poly_n n, -of_poly_sigma f\n");
    printf("  -of_downsample f\n");
}

template<class T>
class PointerList { public:
    vector<T *> L;
    
    ~PointerList() {
        for (int i = 0; i < (int) L.size(); i++) {
            delete L[i];
        }
    }
    
    T *operator[](int i) {
        return L[i];
    }
};

class Pyramid: public PointerList<Array<float> > { public:
    Pyramid(const AnalogyParams &a_params, const Array<float> &a) {
        double max_h = log(a.height());
        double min_h = log(a_params.min_size);
        L.resize(a_params.levels);
        #pragma omp parallel for schedule(dynamic,1)
        for (int i = 0; i < a_params.levels; i++) {
            if (a_params.analogy_verbose) {
                printf("  Create pyramid level %d/%d\n", i, a_params.levels); fflush(stdout);
            }
            double t = i/(a_params.levels-1.0);
            double h = exp(min_h + (max_h-min_h)*t);
            double w = double(a.width())/a.height()*h;
//            double T0 = wall_time();
            int w_current = int(w+0.5);
            int h_current = int(h+0.5);
            if (a_params.analogy_verbose) {
                printf("  Pyramid size: %d x %d, source: %d x %d\n", w_current, h_current, a.width(), a.height()); fflush(stdout);
            }
            L[i] = new Array<float>(imresize(a, w_current, h_current));
            if (a_params.analogy_verbose) {
                printf("  Done creating pyramid level %d/%d\n", i, a_params.levels); fflush(stdout);
            }
//            printf("pyramid level %d (%d x %d): %f\n", i, int(w+0.5), int(h+0.5), wall_time()-T0);
//            Array<float> *current = new Array<float>();
//            imresize_fast(a, *current, int(w+0.5), int(h+0.5));
//            L.push_back(current);
        }
    }
};

double analogy_lerp(double t, double val_at_0, double val_at_1) {
    return val_at_0 + (val_at_1-val_at_0)*t;
}

void get_mean_sigma(const Array<float> &a, double &amean, double &asigma) {
    amean = 0.0;
    asigma = 0.0;
    for (int y = 0; y < a.height(); y++) {
        for (int x = 0; x < a.width(); x++) {
            double c = a(y, x, 0);
            amean += c;
            asigma += c*c;
        }
    }
    int n = a.width() * a.height();
    amean /= n;
    asigma = sqrt((asigma/n) - amean*amean);
}

void lum_remap(Array<float> &a, Array<float> &ap, const Array<float> &b0) {
    Array<float> b(b0);
    
    double amean = 0.0, asigma = 0.0;
    double bmean = 0.0, bsigma = 0.0;
    
    get_mean_sigma(a, amean, asigma);
    get_mean_sigma(b, bmean, bsigma);
    
    for (int y = 0; y < a.height(); y++) {
        for (int x = 0; x < a.width(); x++) {
            a(y, x, 0)  = bsigma / asigma * (a(y, x, 0)  - amean) + bmean;
            ap(y, x, 0) = bsigma / asigma * (ap(y, x, 0) - amean) + bmean;
            
            a(y, x, 0)  = MAX(MIN(a(y, x, 0),  1.0), 0.0);
            ap(y, x, 0) = MAX(MIN(ap(y, x, 0), 1.0), 0.0);
        }
    }
}

void get_composite_descriptor(PatchTableParams *p, AnalogyParams &a_params, const Array<float> &a, const Array<float> &ap, Array<float> &a_desc, Array<float> &ap_desc, Array<float> &desc) {
    if (a_params.analogy_verbose) {
        printf("get_composite_descriptor: a dims: %d, ap dims: %d\n", a.dimensions(), ap.dimensions()); fflush(stdout);
        printf("get_composite_descriptor: gck(a: %dx%dx%d)\n", a.width(), a.height(), a.channels()); fflush(stdout);
    }
    gck<float, float>(a,  a_desc,  p->ndims/2, 0, p->patch_w);
    if (a_params.analogy_verbose) {
        printf("get_composite_descriptor: gck(ap: %dx%dx%d)\n", ap.width(), ap.height(), ap.channels()); fflush(stdout);
    }
    gck<float, float>(ap, ap_desc, p->ndims/2, 0, p->patch_w);
    ASSERT2(a_desc.channels() == p->ndims/2, "expected a_desc channels to be p->ndims/2");
    ASSERT2(ap_desc.channels() == p->ndims/2, "expected ap_desc channels to be p->ndims/2");
    if (a_params.analogy_verbose) {
        printf("get_composite_descriptor: desc, p->ndims=%d\n", p->ndims); fflush(stdout);
    }
    
    desc.resize(a_desc.height(), a_desc.width(), p->ndims);
    for (int y = 0; y < a_desc.height(); y++) {
        for (int x = 0; x < a_desc.width(); x++) {
            for (int c = 0; c < p->ndims; c++) {
                desc(y, x, c) = (c < p->ndims / 2) ? a_desc(y, x, c): (a_params.ap_weight * ap_desc(y, x, c-p->ndims/2));
            }
        }
    }
    if (a_params.analogy_verbose) {
        printf("get_composite_descriptor: done\n"); fflush(stdout);
    }
}

double coherence_weight(const Array<float> &b, const Array<double> &bnn, const Array<float> &a, int ax, int ay, int bx, int by, double upweight_coherence) {

    static double weight_table[9];
    static bool weight_table_init = false;
    if (!weight_table_init) {
        weight_table_init = true;
        for (int i = 0; i < 9; i++) {
            weight_table[i] = 1e-6 + pow(i, upweight_coherence);
        }
    }

    int n = 0;
    int bnn_w = bnn.width();
    int bnn_h = bnn.height();
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx != 0 || dy != 0) {
                int bx_p = bx+dx;
                int by_p = by+dy;
                if (in_bounds(bx_p, bnn_w) && in_bounds(by_p, bnn_h)) {
                    int ax_p = int(bnn(by_p, bx_p, NNF_X));
                    int ay_p = int(bnn(by_p, bx_p, NNF_Y));
                    if (ax_p == ax+dx && ay_p == ay+dy) {
                        n++;
                    }
                }
            }
        }
    }
    
    return weight_table[n];
}

/* Creates image b by using NNF bnn which gives nearest neighbors from b -> a. */
template<int patch_w>
void vote(PatchTableParams *p, AnalogyParams *a_params, Array<float> &b, const Array<double> &bnn, const Array<float> &a, Array<float> &weights) {
    b.resize(bnn.height()+p->patch_w-1, bnn.width()+p->patch_w-1, a.channels());
    if (a_params->vote_overlapping) {
        b.clear(0.0);
    }
    if (a.channels() != 3 || b.channels() != 3) {
        fprintf(stderr, "expected a channels (%d) and b channels (%d) to be 3\n", a.channels(), b.channels());
        exit(1);
    }
    weights.resize(b.height(), b.width(), 1);
    weights.clear(0.0);
    
    #pragma omp parallel for schedule(static, ANALOGY_DEFAULT_PATCH_W)
    for (int by = 0; by < bnn.height(); by++) {
        for (int bx = 0; bx < bnn.width(); bx++) {
            int ax = bnn(by, bx, NNF_X);
            int ay = bnn(by, bx, NNF_Y);
            
            if (!in_bounds(ax, a.width()-p->patch_w+1)) {
#if ANALOGY_ERROR_OUT_OF_BOUNDS
                fprintf(stderr, "expected ax in bounds in nnf, got bx, by=%d, %d, ax, ay=%d, %d, a: %dx%d, patch_w=%d\n", bx, by, ax, ay, a.width(), a.height(), p->patch_w);
                exit(1);
#else
                if (ax < 0) { ax = 0; }
                else if (ax >= a.width()-p->patch_w+1) { ax = a.width()-p->patch_w; }
#endif
            }
            if (!in_bounds(ay, a.height()-p->patch_w+1)) {
#if ANALOGY_ERROR_OUT_OF_BOUNDS
                fprintf(stderr, "expected ay in bounds in nnf, got bx, by=%d, %d, ax, ay=%d, %d, a: %dx%d, patch_w=%d\n", bx, by, ax, ay, a.width(), a.height(), p->patch_w);
                exit(1);
#else
                if (ay < 0) { ay = 0; }
                else if (ay >= a.height()-p->patch_w+1) { ay = a.height()-p->patch_w; }
#endif
            }
            
            if (a_params->vote_overlapping) {
                /* TODO: Could use a falloff such as Gaussian weighting */
                double wpatch = a_params->upweight_coherence ? coherence_weight(b, bnn, a, ax, ay, bx, by, a_params->upweight_coherence): 1.0;
                if (a_params->lum_only) {
#define VOTE_OVERLAPPING(nchannels) \
                    for (int dy = 0; dy < patch_w; dy++) { \
                        float *b_row = b.data + (by+dy) * b.stride[0] + (bx) * b.stride[1]; \
                        float *a_row = a.data + (ay+dy) * a.stride[0] + (ax) * b.stride[1]; \
                        float *weights_row = weights.data + (by+dy) * weights.stride[0] + (bx) * weights.stride[1]; \
                        for (int dx = 0; dx < patch_w; dx++) { \
                            for (int c = 0; c < nchannels; c++) { \
                                b_row[dx*3+c] += a_row[dx*3+c]*wpatch; \
                            } \
                            weights_row[dx] += wpatch; \
                        } \
                    }

//                                b(by+dy, bx+dx, c) += a(ay+dy, ax+dx, c); \
//                            weights(by+dy, bx+dx, 0) += 1; \

                    VOTE_OVERLAPPING(1);
                } else {
                    VOTE_OVERLAPPING(3);
                }
            } else {
                int dy_min = 0, dy_max = p->patch_w;
                int dx_min = 0, dx_max = p->patch_w;
                
                if (bx == 0) { dx_max = p->patch_w/2+1; }
                else if (bx == bnn.width()-1) { dx_min = p->patch_w/2; }
                else { dx_min = p->patch_w/2; dx_max = p->patch_w/2+1; }

                if (by == 0) { dy_max = p->patch_w/2+1; }
                else if (by == bnn.height()-1) { dy_min = p->patch_w/2; }
                else { dy_min = p->patch_w/2; dy_max = p->patch_w/2+1; }
                
#define VOTE_CENTER_PIXEL(nchannels) \
                for (int dy = dy_min; dy < dy_max; dy++) { \
                    for (int dx = dx_min; dx < dx_max; dx++) { \
                        for (int c = 0; c < nchannels; c++) { \
                            b(by+dy, bx+dx, c) = a(ay+dy, ax+dx, c); \
                        } \
                    } \
                }
                if (a_params->lum_only) {
                    VOTE_CENTER_PIXEL(1);
                } else {
                    VOTE_CENTER_PIXEL(3);
                }
            }
        }
    }
    
    if (a_params->vote_overlapping) {
        #pragma omp parallel for schedule(static, ANALOGY_DEFAULT_PATCH_W)
        for (int by = 0; by < b.height(); by++) {
            float *weights_row = weights.data + by * weights.stride[0];
            float *b_row = b.data + by * b.stride[0];
            for (int bx = 0; bx < b.width(); bx++) {
//                double w = weights(by, bx, 0);
                double w = weights_row[bx];
/*
                if (w == 0) {
                    fprintf(stderr, "expected weight nonzero at %d, %d, patch_w=%d, b: %dx%d, bnn: %dx%d, a: %dx%d\n", bx, by, p->patch_w, b.width(), b.height(), bnn.width(), bnn.height(), a.width(), a.height()); exit(1);
                } */
                w = 1.0 / w;
                float *b_pixel = b_row + bx * b.stride[1];

#define VOTE_COMPUTE_MEAN(nchannels) \
                for (int c = 0; c < nchannels; c++) { \
                    b_pixel[c] *= w; \
                }
                if (a_params->lum_only) {
                    VOTE_COMPUTE_MEAN(1);
                } else {
                    VOTE_COMPUTE_MEAN(3);
                }
            }
        }
    }
    if (a_params->lum_only) {
        ASSERT2(b.channels() == 3, "expected 3 channels in lum_only mode");
        if (a_params->log) {
            for (int y = 0; y < b.height(); y++) {
                for (int x = 0; x < b.width(); x++) {
                    b(y, x, 1) = 0;
                    b(y, x, 2) = 0;
                }
            }
        }
    }
}

void log_image(AnalogyParams *a_params, const Array<float> &I0, int level, int iter, string filename, bool is_yuv) {
    static int count = 0;
    if (a_params->log) {
        Array<float> I(I0);
        char full_filename[256];
        sprintf(full_filename, "log_%03d_level%02d_iter%02d_%s", count, level, iter, filename.c_str());
        count++;
        printf("  Before save_color_image\n"); fflush(stdout);
        printf("  I size: %dx%dx%d\n", I.width(), I.height(), I.channels()); fflush(stdout);
        double I_sum = 0.0;
        for (int y = 0; y < I.height(); y++) {
            for (int x = 0; x < I.width(); x++) {
                for (int c = 0; c < I.channels(); c++) {
                    I_sum += I(y, x, c);
                }
            }
        }
        printf("  I sum: %f\n", I_sum); fflush(stdout);
        if (is_yuv) {
            I.yuv2rgb();
        }
        save_color_image<float>(I, full_filename);
        printf("  After save_color_image\n"); fflush(stdout);
    }
}

void list_filenames_matching(AnalogyParams &a_params, const char *filespec_in, const char *filespec_out, vector<string> &L_in, vector<string> &L_out) {
    L_in.clear();
    L_out.clear();
    for (int i = a_params.start_index;; i++) {
        char buf_in[1024];
        char buf_out[1024];
        sprintf(buf_in, filespec_in, i);
        sprintf(buf_out, filespec_out, i-a_params.start_index);
        if (!file_exists(buf_in)) {
            break;
        }
        L_in.push_back(buf_in);
        L_out.push_back(buf_out);
    }
    ASSERT2(L_in.size(), "list_filenames_matching: expected input filename with %d to match at least one file");
}

template<int channels=3>
vector<double> image_histo(const Array<float> &a, int bins_per_channel=8, int step=8) {
    vector<double> ans(bins_per_channel*channels, 0.0);
    for (int y = 0; y < a.height(); y += step) {
        for (int x = 0; x < a.width(); x += step) {
            for (int c = 0; c < channels; c++) {
                double bin = a(y, x, c) * (bins_per_channel - 1e-8);
                int bin_i = int(bin);
                double bin_f = bin-bin_i;
                if (bin_i < 0) {
                    bin_i = 0;
                    bin_f = 0.0;
                } else if (bin_i >= bins_per_channel-1) {
                    bin_i = bins_per_channel-2;
                    bin_f = 1.0;
                }
                ans[bins_per_channel*c+bin_i]   += (1-bin_f);
                ans[bins_per_channel*c+bin_i+1] += bin_f;
            }
        }
    }
    double sum = 0.0;
    for (int i = 0; i < ans.size(); i++) {
        sum += ans[i];
    }
    for (int i = 0; i < ans.size(); i++) {
        ans[i] /= sum;
    }
    
    return ans;
}

class ShotDetector { public:
    double thresh;
    vector<double> prev_histo;
    
    ShotDetector(double thresh_) {
        thresh = thresh_;
    }
    
    bool detect(const Array<float> &frame) {
        if (thresh <= 0) {
            return false;
        }
        bool ans = false;
        
        vector<double> histo = image_histo(frame);
        if (prev_histo.size()) {
            double diff = 0.0;
            for (int i = 0; i < histo.size(); i++) {
                diff += std::abs(histo[i]-prev_histo[i]);
            }
            //printf("ShotDetector diff: %f\n", diff);
            ans = diff > thresh;
        }
        prev_histo = histo;
        return ans;
    }
};

class NNFOpticalFlowCorrector { public:
    AnalogyParams &a_params;
    PatchTableParams *p;
    int a_w, a_h;
    int b_w, b_h;
    cv::Mat flow_mat;
    int frame_index;
    
    NNFOpticalFlowCorrector(AnalogyParams &a_params_, int a_w_, int a_h_, PatchTableParams *p_) :a_params(a_params_) {
        a_w = a_w_;
        a_h = a_h_;
        b_w = b_h = -1;
        frame_index = 0;
        p = p_;
    }
    
    void optical_flow(const Array<float> &b_prev, const Array<float> &b_current) {
        static cv::Ptr<cv::DenseOpticalFlow> dualtv_ptr;
        
        ASSERT2(b_prev.sizes == b_current.sizes, "expected b_prev and b_current sizes to match in optical_flow");
        b_w = b_current.width();
        b_h = b_current.height();
        
        bool use_greyscale = (a_params.of_algo == OPTICAL_FLOW_DUALTV) || (a_params.of_algo == OPTICAL_FLOW_FARNEBACK);
        
        Array<float> b_prev_grey, b_current_grey;
        if (use_greyscale) {
            b_prev_grey.resize(b_prev.height(), b_prev.width(), 1);
            b_current_grey.resize(b_current.height(), b_current.width(), 1);
            for (int y = 0; y < b_prev.height(); y++) {
                for (int x = 0; x < b_prev.width(); x++) {
                    b_prev_grey(y, x, 0) = b_prev(y, x, 0);
                    b_current_grey(y, x, 0) = b_current(y, x, 0);
                }
            }
        }
        const Array<float> &b_prev_ref(use_greyscale ? b_prev_grey: b_prev);
        const Array<float> &b_current_ref(use_greyscale ? b_current_grey: b_current);
        
        cv::Mat b_prev_mat(a_params.of_downsample != 1 ? imdownsample(b_prev_ref, a_params.of_downsample).to_cv8(): b_prev_ref.to_cv8());
        cv::Mat b_current_mat(a_params.of_downsample != 1 ? imdownsample(b_current_ref, a_params.of_downsample).to_cv8(): b_current_ref.to_cv8());
        
        if (a_params.of_algo == OPTICAL_FLOW_SIMPLEFLOW) {
            cv::calcOpticalFlowSF(b_current_mat, b_prev_mat, flow_mat, a_params.of_layers, a_params.of_block_size, a_params.of_max_flow);
        } else if (a_params.of_algo == OPTICAL_FLOW_DUALTV) {
            if (dualtv_ptr.empty()) {
                dualtv_ptr = cv::createOptFlow_DualTVL1();
            }
            
            dualtv_ptr->set("tau", a_params.of_tau);
            dualtv_ptr->set("lambda", a_params.of_lambda);
            dualtv_ptr->set("theta", a_params.of_theta);
            dualtv_ptr->set("nscales", a_params.of_nscales);
            dualtv_ptr->set("epsilon", a_params.of_epsilon);
            dualtv_ptr->set("iterations", a_params.of_iterations > 0 ? a_params.of_iterations: 300);
            
            dualtv_ptr->calc(b_current_mat, b_prev_mat, flow_mat);
        } else if (a_params.of_algo == OPTICAL_FLOW_FARNEBACK) {
            calcOpticalFlowFarneback(b_current_mat, b_prev_mat, flow_mat, a_params.of_pyr_scale, a_params.of_levels, a_params.of_winsize, a_params.of_iterations > 0 ? a_params.of_iterations: 3, a_params.of_poly_n, a_params.of_poly_sigma, 0);
        } else {
            fprintf(stderr, "Unsupported of_algo %d\n", a_params.of_algo); exit(1);
        }
        
        if (a_params.log) {
            char filename[256];
            sprintf(filename, "optical_flow_log%05d.txt", frame_index++);
            Array<float> flow_array(flow_mat);
            save_color_image<float>(flow_array, filename);
        }
    }
    
    void correct(Array<double> &bnn, int a_lo_w, int a_lo_h, int b_lo_w, int b_lo_h, Array<float> &coherence_temporal_sv, int level, Array<float> &accum_dist, bool shot_detected) {
        ASSERT2(b_w > 0, "expected optical_flow() method called before correct()");
        
        if (a_params.analogy_verbose) { printf("in correct\n"); fflush(stdout); }
        if (a_params.analogy_verbose) { printf("bnn sizes: %s\n", vector_to_str_int(bnn.sizes).c_str()); fflush(stdout); }
        Array<double> bnn_copy;
        bnn_copy.assign(bnn);
        coherence_temporal_sv.resize(bnn.height(), bnn.width());
        
        if (a_params.analogy_verbose) { printf("coherence_temporal_sv sizes: %s\n", vector_to_str_int(coherence_temporal_sv.sizes).c_str()); fflush(stdout); }

        double b_to_a_xscale = 1.0*a_w/b_w;
        double b_to_a_yscale = 1.0*a_h/b_h;
        
        double pw2 = ANALOGY_DEFAULT_PATCH_W / 2.0;
        bool is_downsample = a_params.of_downsample != 1.0;
        
        int flow_w = flow_mat.cols, flow_h = flow_mat.rows;
        
        double a_lo_to_hi_xscale = 1.0 * (a_w - 1) / (a_lo_w - 1);
        double a_lo_to_hi_yscale = 1.0 * (a_h - 1) / (a_lo_h - 1);
        
        double b_lo_to_hi_xscale = 1.0 * (b_w - 1) / (b_lo_w - 1);
        double b_lo_to_hi_yscale = 1.0 * (b_h - 1) / (b_lo_h - 1);
        
        double coherence_factor = 1.0/(1+p->coherence_temporal);

        if (a_params.analogy_verbose) { printf("begin loop, flow: %dx%d\n", flow_w, flow_h); fflush(stdout); }
        
        for (int b_lo_y = 0; b_lo_y < bnn.height(); b_lo_y++) {
            for (int b_lo_x = 0; b_lo_x < bnn.width(); b_lo_x++) {
                double b_lo_xcenter = b_lo_x + pw2;
                double b_lo_ycenter = b_lo_y + pw2;
                
                double bx = b_lo_xcenter * b_lo_to_hi_xscale;
                double by = b_lo_ycenter * b_lo_to_hi_yscale;
                
                cv::Vec2f flow_vec;
                if (is_downsample) {
                    double bx_d_real = bx*a_params.of_downsample;
                    double by_d_real = by*a_params.of_downsample;
                    
                    int bx_i = int(bx_d_real);
                    int by_i = int(by_d_real);
                    double bx_f = bx_d_real - bx_i;
                    double by_f = by_d_real - by_i;
                    
                    int bx_i1 = bx_i, bx_i2 = bx_i + 1;
                    int by_i1 = by_i, by_i2 = by_i + 1;

                    if (bx_i1 >= flow_w) { bx_i1 = flow_w-1; }
                    else if (bx_i1 < 0) { bx_i1 = 0; }
                    if (by_i1 >= flow_h) { by_i1 = flow_h-1; }
                    else if (by_i1 < 0) { by_i1 = 0; }

                    if (bx_i2 >= flow_w) { bx_i2 = flow_w-1; }
                    else if (bx_i2 < 0) { bx_i2 = 0; }
                    if (by_i2 >= flow_h) { by_i2 = flow_h-1; }
                    else if (by_i2 < 0) { by_i2 = 0; }
                    
                    cv::Vec2f vec11(flow_mat.at<cv::Vec2f>(by_i1, bx_i1));
                    cv::Vec2f vec12(flow_mat.at<cv::Vec2f>(by_i1, bx_i2));
                    cv::Vec2f vec21(flow_mat.at<cv::Vec2f>(by_i2, bx_i1));
                    cv::Vec2f vec22(flow_mat.at<cv::Vec2f>(by_i2, bx_i2));
                    
                    cv::Vec2f vec1 = vec11 + (vec12-vec11) * bx_f;
                    cv::Vec2f vec2 = vec21 + (vec22-vec21) * bx_f;
                    
                    flow_vec = vec1 + (vec2-vec1) * by_f;
                    flow_vec = flow_vec / a_params.of_downsample;
                } else {
                    int bx_d = int(bx+0.5);
                    int by_d = int(by+0.5);
                    if (bx_d < 0) { bx_d = 0; }
                    else if (bx_d >= flow_w) { bx_d = flow_w-1; }
                    if (by_d < 0) { by_d = 0; }
                    else if (by_d >= flow_h) { by_d = flow_h-1; }
                    
                    flow_vec = flow_mat.at<cv::Vec2f>(by_d, bx_d);
                }
                float flow_x = flow_vec[0];
                float flow_y = flow_vec[1];
                
#if ANALOGY_MAX_FLOW_DIST
                if (level == a_params.levels - 1) {
                    float flow_dist = std::abs(flow_x) + std::abs(flow_y);
                    accum_dist(b_lo_y, b_lo_x) += flow_dist;
                }
#endif
                
                double b_lo_x_p = b_lo_x + flow_x / b_lo_to_hi_xscale;
                double b_lo_y_p = b_lo_y + flow_y / b_lo_to_hi_yscale;
                
                bool use_coh = true;
                if (shot_detected) { use_coh = false; }
                
                int b_lo_x_p_i = int(b_lo_x_p+0.5);
                int b_lo_y_p_i = int(b_lo_y_p+0.5);
                if (b_lo_x_p_i < 0) { b_lo_x_p_i = 0; use_coh = false; }
                else if (b_lo_x_p_i >= bnn.width()) { b_lo_x_p_i = bnn.width()-1; use_coh = false; }
                if (b_lo_y_p_i < 0) { b_lo_y_p_i = 0; use_coh = false; }
                else if (b_lo_y_p_i >= bnn.height()) { b_lo_y_p_i = bnn.height()-1; use_coh = false; }

#if ANALOGY_MAX_FLOW_DIST
                int b_hi_x = (b_lo_xcenter * b_lo_to_hi_xscale) - pw2;
                int b_hi_y = (b_lo_ycenter * b_lo_to_hi_yscale) - pw2;
                if (b_hi_x < 0) { b_hi_x = 0; }
                else if (b_hi_x >= accum_dist.width()) { b_hi_x = accum_dist.width()-1; }
                if (b_hi_y < 0) { b_hi_y = 0; }
                else if (b_hi_y >= accum_dist.height()) { b_hi_y = accum_dist.height()-1; }

                float flow_dist_t = accum_dist(b_hi_y, b_hi_x) / a_params.max_flow_dist;
                if (flow_dist_t > 1.0 ) { flow_dist_t = 1.0; }
                
                double w_t = 1.0/(1+analogy_lerp(flow_dist_t, p->coherence_temporal, 0.0));
                if (!use_coh) { w_t = 1.0; }
#else
                double w_t = coherence_factor;
                if (!use_coh) { w_t = 1.0; }
#endif
                if (a_params.lock && level < a_params.lock) {
                    if (use_coh) {
                        w_t = 0.0;      /* Skip a patch search at this level */
                    }
                }
                coherence_temporal_sv(b_lo_y, b_lo_x) = w_t;
                
                int ax = bnn_copy(b_lo_y_p_i, b_lo_x_p_i, NNF_X);
                int ay = bnn_copy(b_lo_y_p_i, b_lo_x_p_i, NNF_Y);
                
                bnn(b_lo_y, b_lo_x, NNF_X) = ax;
                bnn(b_lo_y, b_lo_x, NNF_Y) = ay;
            }
        }
    }
};

void analogy(int argc, char **argv) {
    ASSERT2(argc >= 4, "analogy expected 4 arguments");
    
    const char *filename_in = argv[2];
    const char *filename_out = argv[3];
    vector<string> L_in({string(filename_in)});
    vector<string> L_out({string(filename_out)});
    
    double T0_analogy = wall_time();
    
    double T_synthesis = 0.0;
    double T_io = 0.0;
    double T_table = 0.0;
    double T_search = 0.0;
    double T_get_composite_descriptor = 0.0;
    double T_vote = 0.0;
    double T_upsample = 0.0;
    double T_pyr = 0.0;
    double T_optical_flow = 0.0;
    double T_shot_detect = 0.0;
    
    AnalogyParams a_params(argc, argv);
    
    if (strstr(filename_in, "%")) {
        list_filenames_matching(a_params, filename_in, filename_out, L_in, L_out);
    }
    
    Array<float> a(load_color_image<float>(argv[0]));
    Array<float> ap(load_color_image<float>(argv[1]));
    Array<float> b(load_color_image<float>(L_in[0].c_str()));

    double T_begin_table = wall_time();
    T_io += T_begin_table - T0_analogy;
    
    a.rgb2yuv();
    ap.rgb2yuv();
    b.rgb2yuv();
    
    lum_remap(a, ap, b);

    log_image(&a_params, a, 0, 0, "a_remap.png", true);
    log_image(&a_params, ap, 0, 0, "ap_remap.png", true);

    PatchTableParams *p = new PatchTableParams();
    p->verbose = 0;
    p->ndims = 8;
    p->coherence_spatial = 80;
    p->coherence_temporal = 1.0;
    p->convert_colorspace = false;
    p->descriptor_padded = true;
    p->set_from_switches(argc, argv);
    ASSERT2(p->ndims % 2 == 0, "analogy: expected ndims to be even");
    p->nchroma = 0;
    
    cv::setNumThreads(p->threads);
    
    PatchTableParams p0(*p);
    p0.ndims /= 2;
    
    Pyramid a_pyr(a_params, a);
    Pyramid ap_pyr(a_params, ap);
    
    PointerList<PatchTable<> > a_table;
    PatchTable<> *a0_table = NULL;          /* Initial lookup table, with half as many dimensions, because bp is not yet known. */

    /* Build table at each pyramid level (a_table) */
    
    for (int level = 0; level < a_params.levels; level++) {
        if (a_params.analogy_verbose) {
            printf("Build table, level %d: a: %dx%d, ap: %dx%d\n", level, a_pyr[level]->width(), a_pyr[level]->height(),
                                                                          ap_pyr[level]->width(), ap_pyr[level]->height()); fflush(stdout);
        }
        Array<float> a_desc, ap_desc, desc;
        if (a_params.analogy_verbose) { printf("Before get_composite_descriptor\n"); fflush(stdout); }
        get_composite_descriptor(p, a_params, *a_pyr[level], *ap_pyr[level], a_desc, ap_desc, desc);
        
        if (a_params.analogy_verbose) { printf("Before a_table.push_back(PatchTable)\n"); fflush(stdout); }
        p->is_descriptor = true;
        a_table.L.push_back(new PatchTable<float, float>(p, desc, (Array<int> *) NULL));
        if (a_params.analogy_verbose) { printf("After a_table.push_back(PatchTable)\n"); fflush(stdout); }
        
        if (level == 0) {
            p->is_descriptor = false;
            a0_table = new PatchTable<float, float>(&p0, *a_pyr[level], (Array<int> *) NULL);
            if (a_params.analogy_verbose) { printf("After a0_table.push_back(PatchTable)\n"); fflush(stdout); }
        }
        
        if (a_params.log) { printf("Saving a.png at current level\n"); fflush(stdout); }
        log_image(&a_params, *a_pyr[level],  level, 0, "a.png", true);
        if (a_params.log) { printf("Saving ap.png at current level\n"); fflush(stdout); }
        log_image(&a_params, *ap_pyr[level], level, 0, "ap.png", true);
        if (a_params.log) { printf("Saved all a images at current level\n"); fflush(stdout); }
    }
    
    printf("Synthesis\n"); fflush(stdout);
    
    T_table += wall_time() - T_begin_table;
    
    /* Synthesis */

    Array<double> b_prev;
    ShotDetector shot_detector(a_params.shot_detect);
    
    Array<float> accum_dist;             /* Accumulated distance of warping due to optical flow */
#if ANALOGY_MAX_FLOW_DIST
    accum_dist.resize(b.height() - p->patch_w+1, b.width() - p->patch_w+1);
    accum_dist.clear(0.0);
#endif
    
    vector<Array<double> > bnn_pyr;       /* Maps b patch coords to a patch coords (copies pixels from ap into b). */
    for (int level = 0; level < a_params.levels; level++) {
        bnn_pyr.push_back(Array<double>(1, 1));
    }
    
    NNFOpticalFlowCorrector *corrector = a_params.optical_flow ? new NNFOpticalFlowCorrector(a_params, a.width(), a.height(), p): NULL;
    
    double coherence_spatial0 = p->coherence_spatial;
    double upweight_coherence0 = a_params.upweight_coherence;
    vector<int> shot_detected_frames;
    
    for (int iframe = 0; iframe < (int) L_out.size(); iframe++) {
        printf("\n");
        int nframes = iframe ? iframe: 1;
        printf("----------------------------------------------------------------------\n");
        printf("Synthesizing frame %d/%d (%.2f%%), synthesis time/frame: %f, search time/frame: %f\n", iframe+1, (int) L_out.size(), iframe*100.0/(L_out.size()), T_synthesis/nframes, T_search/nframes);
        printf("----------------------------------------------------------------------\n");
        if (a_params.analogy_verbose) { printf("\n"); }
        
        if (a_params.optical_flow) {
            b_prev.assign(b);
        }
        bool shot_detected = false;
        if (iframe) {
            double T0_io = wall_time();
            b.assign(load_color_image<float>(L_in[iframe].c_str()));
            T_io += wall_time() - T0_io;
        }

        double T0_synthesis = wall_time();
        
        if (iframe) {
            double T0_shot_detect = wall_time();
            shot_detected = shot_detector.detect(b);
            T_shot_detect += wall_time() - T0_shot_detect;
            b.rgb2yuv();
        }
        if (a_params.analogy_verbose) {
            if (shot_detected) {
                shot_detected_frames.push_back(iframe);
            }
            printf("  shot_detected_frames=%s\n", vector_to_str_int(shot_detected_frames).c_str());
        }
        
        if (a_params.optical_flow) {
            double T0_optical_flow = wall_time();
            corrector->optical_flow(b_prev, b);
            T_optical_flow += wall_time()-T0_optical_flow;
        }
        
        double T0_pyr = wall_time();
        Pyramid b_pyr(a_params, b);
        if (a_params.analogy_verbose) {
            printf("B pyramid sizes:\n");
            for (int level = 0; level < a_params.levels; level++) {
                printf("  level %d: %dx%d\n", level, b_pyr[level]->width(), b_pyr[level]->height());
            }
            printf("\n");
        }
        T_pyr += wall_time()-T0_pyr;

        Array<float> bp;
        Array<float> b_desc, bp_desc, desc;
        Array<float> weights;
        Array<float> coherence_temporal_sv;
        coherence_temporal_sv.resize(1, 1);
        bool coherence_temporal_init = false;
        
        bp.resize(1, 1, 1);
        
        for (int level = 0; level < a_params.levels; level++) {
            p->coherence_spatial = coherence_spatial0;
            a_params.upweight_coherence = upweight_coherence0;
            
            if (level >= a_params.levels - a_params.fine_levels) {
                if (a_params.fine_coherence_spatial >= 0) {
                    p->coherence_spatial = a_params.fine_coherence_spatial;
                }
                if (a_params.fine_upweight_coherence >= 0) {
                    a_params.upweight_coherence = a_params.fine_upweight_coherence;
                }
            }

            if (a_params.analogy_verbose) {
                printf("\nSynthesis, level %d, begin\n", level);
            }
        
            Array<double> bnn_prev;
            Array<double> &bnn(bnn_pyr[level]);
/*
            if (bnn.width() == 1 && bnn.height() == 1) {
                bnn.resize(b_pyr[level]->height()-p->patch_w+1, b_pyr[level]->width()-p->patch_w+1, 3);
                bnn.clear(0.0);
            }
 */
            if (a_params.optical_flow) {
                double T0_optical_flow = wall_time();
                if (a_params.analogy_verbose) {
                    printf("\nSynthesis, level %d, before correct\n", level);
                }
                corrector->correct(bnn, ap_pyr[level]->width(), ap_pyr[level]->height(), b_pyr[level]->width(), b_pyr[level]->height(), coherence_temporal_sv, level, accum_dist, shot_detected);
                if (a_params.analogy_verbose) {
                    printf("\nSynthesis, level %d, after correct, coherence_temporal_init=%d, coherence_temporal_sv sizes: %s\n", level, int(coherence_temporal_init), vector_to_str_int(coherence_temporal_sv.sizes).c_str());
                }
                T_optical_flow += wall_time()-T0_optical_flow;
            }
            bnn_prev.assign(bnn);
            
            if (a_params.analogy_verbose) {
                printf("\nSynthesis, level %d, ap: %dx%dx%d\n", level, ap_pyr[level]->width(), ap_pyr[level]->height(), ap_pyr[level]->channels());
            }

            if (a_params.log) { printf("Saving b.png at current level\n"); fflush(stdout); }
            log_image(&a_params, *b_pyr[level],  level, 0, "b.png", true);

            for (int iter = 0; iter < a_params.correction_iters; iter++) {
                if (a_params.analogy_verbose) {
                    printf("  Iteration %d, ap: %dx%dx%d, bp: %dx%dx%d\n", iter, ap_pyr[level]->width(), ap_pyr[level]->height(), ap_pyr[level]->channels(), bp.width(), bp.height(), bp.channels());
                }
                
                /* Update NNF bnn */
                Array<float> *coherence_temporal_sv_p = NULL;
                if (coherence_temporal_init && coherence_temporal_sv.width() > 1 && coherence_temporal_sv.height() > 1) {
                    coherence_temporal_sv_p = &coherence_temporal_sv;
                }
                if (iframe % a_params.update_every == 0) {
                    if (level == 0 && iter == 0) {
                        if (a_params.analogy_verbose) { printf("  Doing initial lookup, coherence_temporal_sv_p=%p, coherence_temporal_sv sizes=%s\n", coherence_temporal_sv_p, coherence_temporal_sv_p ? vector_to_str_int(coherence_temporal_sv.sizes).c_str(): ""); }
                        p->is_descriptor = false;
                        double T0_search = wall_time();
                        a0_table->lookup(*b_pyr[0], bnn, bnn_prev.dimensions() == 3 ? &bnn_prev: NULL, coherence_temporal_sv_p);
                        T_search += wall_time()-T0_search;
                        if (a_params.analogy_verbose) { printf("  bnn: %dx%dx%d\n", bnn.width(), bnn.height(), bnn.channels()); }
                    } else {
                        if (a_params.analogy_verbose) { printf("  Doing iterative lookup -- computing descriptor\n"); }
                        double T0_search = wall_time();
                        get_composite_descriptor(p, a_params, *b_pyr[level], bp, b_desc, bp_desc, desc);
                        if (a_params.analogy_verbose) {
                            printf("  Doing iterative lookup -- lookup(), bp: %dx%dx%d, desc: %dx%dx%d, p=%p, patch_w=%d, use_coherence_temporal_p=%p, coherence_temporal_sv sizes=%s\n", bp.width(), bp.height(), b.channels(), desc.width(), desc.height(), desc.channels(), (void *) p, p->patch_w, coherence_temporal_sv_p, coherence_temporal_sv_p ? vector_to_str_int(coherence_temporal_sv.sizes).c_str(): "");
                        }
                        T_get_composite_descriptor += wall_time()-T0_search;
                        
                        p->is_descriptor = true;
                        a_table[level]->lookup(desc, bnn, bnn_prev.dimensions() == 3 ? &bnn_prev: NULL, coherence_temporal_sv_p);

#if ANALOGY_MAX_FLOW_DIST
                        if (a_params.analogy_verbose) { printf("  done lookup, setting accum_dist\n"); }
                        for (int y = 0; y < bnn.height(); y++) {
                            for (int x = 0; x < bnn.width(); x++) {
                                if (int(bnn(y, x, NNF_X)) != int(bnn_prev(y, x, NNF_X)) ||
                                    int(bnn(y, x, NNF_Y)) != int(bnn_prev(y, x, NNF_Y))) {
                                    accum_dist(y, x) = 0;
                                }
                            }
                        }
#endif
                        T_search += wall_time()-T0_search;
                        if (a_params.analogy_verbose) { printf("  bnn: %dx%dx%d\n", bnn.width(), bnn.height(), bnn.channels()); }
                    }
                }
                
                /* Update bp image */
                if (a_params.analogy_verbose) {
                    printf("  Before vote, bp: %dx%dx%d, bnn: %dx%dx%d, ap: %dx%dx%d\n", bp.width(), bp.height(), bp.channels(), bnn.width(), bnn.height(), bnn.channels(), ap_pyr[level]->width(), ap_pyr[level]->height(), ap_pyr[level]->channels());
                }
                double T0_vote = wall_time();
                vote<ANALOGY_DEFAULT_PATCH_W>(p, &a_params, bp, bnn, *ap_pyr[level], weights);
                T_vote += wall_time()-T0_vote;
                if (a_params.analogy_verbose) { printf("  After vote, bp: %dx%dx%d\n", bp.width(), bp.height(), bp.channels()); }

                log_image(&a_params, bnn, level, iter, "bnn.pfm", false);
                if (a_params.log) { printf("  Before log bp.pfm\n"); fflush(stdout); }
                if (a_params.log) {
                    Array<float> bresized(imresize(b, bp.width(), bp.height()));
                    if (a_params.lum_only) {
                        /* Put color channels back */
                        ASSERT2(bresized.sizes == bp.sizes, "expected bresized and bp to be same size");
                        ASSERT2(bresized.channels() == 3, "expected bresized channels to be 3");
                        for (int y = 0; y < bp.height(); y++) {
                            for (int x = 0; x < bp.width(); x++) {
                                bp(y, x, 1) = bresized(y, x, 1);
                                bp(y, x, 2) = bresized(y, x, 2);
                            }
                        }
                    }
                }
                log_image(&a_params, bp, level, 0, "bp.pfm", false);
                if (a_params.log) { printf("  After log bp.pfm, saving bp.png\n"); fflush(stdout); }
                log_image(&a_params, bp, level, 0, "bp.png", true);
                if (a_params.log) { printf("  After log bp.png\n"); fflush(stdout); }

                if (iter+1 < a_params.correction_iters) {
                    bnn_prev.assign(bnn);
                }
            }
            
            /* Upsample bp image */
            if (level+1 < a_params.levels) {
                /* TODO: Upsample using NNF */
                double T0_upsample = wall_time();
                if (a_params.analogy_verbose) {
                    printf("  Before resize bp: %dx%dx%d\n", bp.width(), bp.height(), bp.channels()); fflush(stdout);
                }
                bp = imresize(bp, b_pyr[level+1]->width(), b_pyr[level+1]->height());  /* TODO: This could be done in-place for speed */
                if (a_params.analogy_verbose) {
                    printf("  After resize bp\n"); fflush(stdout);
                }
//                imresize_fast(bp, bp, b_pyr[level+1]->width(), b_pyr[level+1]->height());  /* TODO: This could be done in-place for speed */
                T_upsample += wall_time()-T0_upsample;
            }
            coherence_temporal_init = true;
        }

        if (a_params.lum_only) {
            /* Put color channels back */
            ASSERT2(bp.sizes == b.sizes, "expected b and bp to be same size");
            ASSERT2(bp.channels() == 3, "expected bp channels to be 3");
            for (int y = 0; y < bp.height(); y++) {
                for (int x = 0; x < bp.width(); x++) {
                    bp(y, x, 1) = b(y, x, 1);
                    bp(y, x, 2) = b(y, x, 2);
                }
            }
        }
        bp.yuv2rgb();
        
        double T_end_synthesis = wall_time();
        T_synthesis += T_end_synthesis - T0_synthesis;
        save_color_image<float>(bp, L_out[iframe].c_str());
        T_io += wall_time() - T_end_synthesis;
        
        double T_end = wall_time();
        
        printf("  Analogy time:                   %f\n", T_end - T0_analogy);
        printf("    Build table:                  %f\n", T_table);
        printf("    Synthesis:                    %f\n", T_synthesis);
        printf("      Shot detection:             %f (%.1f%%)\n", T_shot_detect, T_shot_detect*100.0/T_synthesis);
        printf("      Pyramid:                    %f (%.1f%%)\n", T_pyr, T_pyr*100.0/T_synthesis);
        printf("      Search:                     %f (%.1f%%)\n", T_search, T_search*100.0/T_synthesis);
        printf("        get_composite_descriptor: %f (%.1f%%)\n", T_get_composite_descriptor, T_get_composite_descriptor*100.0/T_synthesis);
        printf("      Vote:                       %f (%.1f%%)\n", T_vote, T_vote*100.0/T_synthesis);
        printf("      Upsample:                   %f (%.1f%%)\n", T_upsample, T_upsample*100.0/T_synthesis);
        printf("      Optical flow:               %f (%.1f%%)\n", T_optical_flow, T_optical_flow*100.0/T_synthesis);
        printf("    Disk IO:                      %f\n", T_io);
    }

    delete p;
    delete a0_table;
}

