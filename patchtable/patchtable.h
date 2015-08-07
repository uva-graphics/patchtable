
/* ------------------------------------------------------------------
   PatchTable header. See patchtable_main.cpp for example usage
   ------------------------------------------------------------------ */
#ifndef _patchtable_h
#define _patchtable_h

#include <algorithm>
#include <memory>
#include <boost/pending/disjoint_sets.hpp>
#include <unordered_map>
using std::sort;
using std::shared_ptr;
using std::unordered_map;
using std::make_heap;
using std::push_heap;
using std::pop_heap;

#define TABLE_DEBUG                   0        /* Extra debug checks and prints (should be turned off for checked in code) */
#define TABLE_DEBUG_SERIALIZE         0        /* Same for serialization */

#define TABLE_VERBOSE                 0        /* Extra verbosity (should be turned off for checked in code) */
#define TABLE_EXTRA_VERBOSE           0        /* Even more verbosity (should be turned off for checked in code) */
#define PM_VERBOSE                    0        /* Verbose in PatchMatch */

#define TABLE_OPENMP                  1        /* Use OpenMP (multithreading) */

#define TABLE_SSE                     1        /* Use streaming SIMD instructions (for x86 architecture), a slight optimization. */
#define TABLE_INT                     0        /* Use integer mode for table (does not work) */

#define TABLE_OPTIMIZE_DT             1        /* Optimize dt calculation in exchange for using more memory. Must be 0 if TABLE_INFINITY_NORM is 1 (unless TABLE_DT_REGULAR is 1, in which case this option is ignored). */
#define TABLE_RANDOMIZE_ADD           0        /* Randomize when adding patches in table */
#define TABLE_OPTIMIZE_REGULAR        0        /* Always uses regular_grid mode and also optimizes dt computation */
#define TABLE_DT_REGULAR              1        /* Run dt assuming regular grid (even if grid is not in fact regular) */
#define TABLE_INFINITY_NORM           1        /* Use infinity norm for table (should be 1 for best results) */

#define TABLE_POPULATE_FAST           1        /* Whether to populate fast. Should be turned off for lower error. */

//#define TABLE_RANDOMIZE_LOOKUP        0        /* Randomize lookup tables */

#define TABLE_TREECANN_ANN            1        /* Use ANN search library for TreeCANN implementation */

#define TABLE_ENABLE_FLANN            1        /* Enable FLANN search library */
#define TABLE_ENABLE_ANN              1        /* Enable ANN search library (this should be 1 if TABLE_TREECANN_ANN is 1) */

#define TABLE_COUNT_BINS              0        /* Should be 0 for checked in code */

#define TABLE_PROP_FAST               1        /* Use fast summed-area table propagation (should be 1 for checked in code) */

#define TABLE_ALLOW_COMPARE           0
#define TABLE_PROP_EXACT              0        /* Use exact distance in table lookup */
#define TABLE_PATCHMATCH_EXACT        0        /* Use exact distance in PatchMatch precomputation. Should be 0 (gives lower error). */

#define TABLE_DT_DOWNSAMPLE           1        /* Enable downsample dt */

#define TABLE_PATCH_W                 8        /* TODO This should actually be templated */
#define TABLE_IMAGE_CHANNELS          3

#define TABLE_ALLOW_QUERY_STEP        1

#define TABLE_SIZE_STDDEV             0        /* Make count be proportional to standard deviation along each dimension */

#define TABLE_PROFILE                 0        /* Additional profiling of time */

#define TABLE_KCOHERENCE_STEP         1        /* Enable kcoherence_step mode */

#define TABLE_MULTI_TABLES            1        /* Enable ntables parameter */

#define TABLE_KCOHERENCE_ENRICH       0        /* Enable kcoherence_enrich mode */
#define TABLE_SAVE_KCOHERENCE         0        /* Enable save_kcoherence mode */

#define TABLE_KCOHERENCE_TRIANGLE     1        /* Accelerate kcoherence with triangle inequality */

#define TABLE_PRODUCT_QUANTIZE_ACTUAL_DIST 1   /* Use actual distance in distance transform */

#if TABLE_ENABLE_FLANN
#define INF 1E20
//#define BUILD_UNIX 1
//#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include <flann/flann.hpp>
#endif

#if BUILD_DEBUG
#undef TABLE_DEBUG
#undef TABLE_VERBOSE
#undef TABLE_EXTRA_VERBOSE
#define TABLE_DEBUG 1
#define TABLE_VERBOSE 0
#define TABLE_EXTRA_VERBOSE 0
#endif

#include "patchtable_params.h"
#include "gck.h"
#include "patch_pca.h"

#if TABLE_ENABLE_ANN
#include <ANN/ANN.h>
#endif

#if TABLE_SSE
#include <mmintrin.h>
#endif

#if TABLE_OPENMP
#include <omp.h>
#endif

#define TABLE_MASK_SHIFT 30
#define TABLE_LO_MASK ((1<<TABLE_MASK_SHIFT)-1)
#define TABLE_HI_MASK (1<<TABLE_MASK_SHIFT)
#define TABLE_UNUSED (TABLE_LO_MASK)

#if TABLE_DT_DOWNSAMPLE
#include "dt_downsample.h"
#endif

#define TABLE_DEFAULT_ITYPE int32_t

#define TABLE_FOR_ALLOWED_PATCHES_ARRAY_STEP(wh0, step) \
for (int y = gck_ymin(wh0); y < gck_ymax(wh0); y += step) { \
    for (int x = gck_xmin(wh0); x < gck_xmax(wh0); x += step) { \
        if (allowed_patches) { \
            int xpatch = x-gck_xmin(wh0); \
            int ypatch = y-gck_ymin(wh0); \
            if ((*allowed_patches)(ypatch, xpatch) != p->allowed_index) { \
                continue; \
            } \
        }

#define TABLE_FOR_ALLOWED_PATCHES_ARRAY(wh0) \
    TABLE_FOR_ALLOWED_PATCHES_ARRAY_STEP(wh0, 1)

#define TABLE_FOR_ALLOWED_PATCHES_STEP(step) \
    TABLE_FOR_ALLOWED_PATCHES_ARRAY_STEP(wh0, step)

#define TABLE_FOR_ALLOWED_PATCHES() TABLE_FOR_ALLOWED_PATCHES_STEP(1)

#if TABLE_CLUSTER_KMEANS
#define TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER() TABLE_FOR_ALLOWED_PATCHES_ARRAY(choose_wh0)
#else
#define TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER() TABLE_FOR_ALLOWED_PATCHES()
#endif

#define TABLE_END_FOR_ALLOWED_PATCHES() } }

#if TABLE_PRODUCT_QUANTIZE
#include "product_quantize.h"
#endif

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,    0> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,    0> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<int,    1> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<int,    1> >,
  boost::mpl::pair<uint32_t, boost::mpl::integral_c<int,    1> >,
  boost::mpl::pair<int8_t,   boost::mpl::integral_c<int,    1> >,
  boost::mpl::pair<int16_t,  boost::mpl::integral_c<int,    1> >,
  boost::mpl::pair<int32_t,  boost::mpl::integral_c<int,    1> > > epsilonByType;

#define EPSILON_BY_TYPE(T) (boost::mpl::at<epsilonByType, T>::type::value == 0 ? 1e-4: boost::mpl::at<epsilonByType, T>::type::value)

#define XY_TO_INT(x, y) ((x)|((y)<<16))
#define INT_TO_X(v) ((v)&(65535))
#define INT_TO_Y(v) ((v)>>16)

#if TABLE_CLUSTER_KMEANS
#define XY_TO_INT_OPTIONAL_CLUSTER(x, y) (p->cluster_kmeans ? (x): XY_TO_INT(x, y))
#else
#define XY_TO_INT_OPTIONAL_CLUSTER(x, y) XY_TO_INT(x, y)
#endif

#define DIST_INFINITY 100000

#define NNF_X    0
#define NNF_Y    1
#define NNF_DIST 2

#if TABLE_OPTIMIZE_REGULAR
#define dist_t int
#else
#define dist_t float
#endif

#define TABLE_TRY_DT() \
                            int vsrc = table(src_index); \
                            if (vsrc != TABLE_UNUSED) { \
                                int xsrc = INT_TO_X(vsrc); \
                                int ysrc = INT_TO_Y(vsrc); \
                                dist_t dcurrent = part->template patch_dist_to_grid<1>(wh0, xsrc, ysrc, grid_index, dbest); \
                                if (dcurrent < dbest) { \
                                    dbest = dcurrent; \
                                    xbest = xsrc; \
                                    ybest = ysrc; \
                                } \
                            }

#define DT_MODE_MANHATTAN       0
#define DT_MODE_EXPTIME         1
#define DT_MODE_SEMI_EUCLIDEAN  2

#define DT_ALGO_RASTER            0
#define DT_ALGO_PROP              1
#define DT_ALGO_BRUTE             2
#define DT_ALGO_KDTREE            3
#define DT_ALGO_DOWNSAMPLE        4
#define DT_ALGO_HYBRID            5
#define DT_ALGO_EUCLIDEAN         6

#define KCOHERENCE_ALGO_FLANN     0
#define KCOHERENCE_ALGO_PM        1
#define KCOHERENCE_ALGO_ANN       2
#define KCOHERENCE_ALGO_TREECANN  3

//#define PATCH_PARTITION_NMAP 1024
#define PATCH_PARTITION_NMAP 16384

#define TABLE_DEFAULT_MIN_DIST  2

/* ------------------------------------------------------------------------------------
   Parameters, which can be initialized from command-line switches
   ------------------------------------------------------------------------------------ */

map<string, string> parse_switches(int argc, char **argv);

#define PATCHTABLE_COLORSPACE_YUV 0
#define PATCHTABLE_COLORSPACE_LAB 1

#define LOOKUP_ALGO_TABLE      0
#define LOOKUP_ALGO_BRUTE      1
#define LOOKUP_ALGO_KDTREE     2
#define LOOKUP_ALGO_PM         3
#define LOOKUP_ALGO_TREECANN   4

#define TABLE_DIM_ALGO_WH      0
#define TABLE_DIM_ALGO_PCA     1

void print_switches();

/* ------------------------------------------------------------------------------------------------
   Serialization of basic data-types
   ------------------------------------------------------------------------------------------------ */

#define LOAD_SAVE_DTYPE_HEADER(dtype) \
void save(FILE *f, dtype v); \
void load(FILE *f, dtype &v); \
void save(FILE *f, const vector<dtype> &v); \
void load(FILE *f, vector<dtype> &v); \
void save(FILE *f, const Array<dtype> &v); \
void load(FILE *f, Array<dtype> &v);

LOAD_SAVE_DTYPE_HEADER(int32_t);
LOAD_SAVE_DTYPE_HEADER(char);
LOAD_SAVE_DTYPE_HEADER(int8_t);
LOAD_SAVE_DTYPE_HEADER(int16_t);
LOAD_SAVE_DTYPE_HEADER(int64_t);
LOAD_SAVE_DTYPE_HEADER(float);
LOAD_SAVE_DTYPE_HEADER(double);

void save(FILE *f, const string &s);
void save(FILE *f, const char *s);
void load(FILE *f, string &s);
void load_id(FILE *f, const string &s);

/* ------------------------------------------------------------------------------------------------
   Basic matching algorithms
   ------------------------------------------------------------------------------------------------ */

#define PATCH_DIST_CALC_BVEC0(wh0) \
        const real *bvec0 = &wh0.get_nearest((p->patch_w-1), (p->patch_w-1), 0);

#define PATCH_DIST_CALC_AVEC_NO_DECL_COORD(wh, x, y) \
        int ycurrent = y+p->patch_w-1; \
        int xcurrent = x+p->patch_w-1; \
        avec = &wh.get_nearest(ycurrent, xcurrent, 0);

#define PATCH_DIST_CALC_AVEC_NO_DECL(wh) \
        PATCH_DIST_CALC_AVEC_NO_DECL_COORD(wh, x, y);

#define PATCH_DIST_CALC_AVEC_COORD(wh, x, y) \
        const real *avec; \
        PATCH_DIST_CALC_AVEC_NO_DECL_COORD(wh, x, y);

#define PATCH_DIST_CALC_AVEC(wh) \
        const real *avec; \
        PATCH_DIST_CALC_AVEC_NO_DECL(wh);

#define PATCH_DIST_CALC_EXACT_AVEC_NO_DECL_COORD(a, x, y) \
        avec = &a.get_nearest(y, x, 0);

#define PATCH_DIST_CALC_EXACT_AVEC_NO_DECL(a) \
        PATCH_DIST_CALC_EXACT_AVEC_NO_DECL_COORD(a, x, y)

#define PATCH_DIST_CALC_EXACT_AVEC(a) \
        const real *avec; \
        PATCH_DIST_CALC_EXACT_AVEC_NO_DECL(a);

#define PATCH_DIST_CALC_EXACT_AVEC_COORD(a, x, y) \
        const real *avec; \
        PATCH_DIST_CALC_EXACT_AVEC_NO_DECL_COORD(a, x, y);

#define PATCH_DIST_CALC_AVEC_EXACT_OR_APPROX(wh_a, a, x, y, is_descriptor) \
        const real *avec; \
        if (is_descriptor) { \
            PATCH_DIST_CALC_AVEC_NO_DECL_COORD(wh_a, x, y); \
        } else { \
            PATCH_DIST_CALC_EXACT_AVEC_NO_DECL_COORD(a, x, y); \
        }

#define TABLE_PREPARE_B_DIST() \
    const real *avec = p->is_descriptor ? &wh0.get_nearest(y, x, 0): &b0.get_nearest(y-gck_ymin(wh0), x-gck_xmin(wh0), 0); \
    int is_descriptor = p->is_descriptor; \
    const Array<real> &a(b0);

#define TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, is_descriptor) \
        (is_descriptor ? patch_dist_approx<real, TABLE_NCHANNELS>(wh0, avec, bvec0, xsrc, ysrc): \
                         patch_dist_exact<real>(a, avec, b0, xsrc, ysrc))

#define TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc) \
        TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, is_descriptor)

#define TABLE_PRECALC_A_PATCH() \
        const real *avec; \
        if (is_descriptor) { \
            PATCH_DIST_CALC_AVEC_NO_DECL(a_wh); \
        } else { \
            PATCH_DIST_CALC_EXACT_AVEC_NO_DECL(a); \
        }

template<class real, int patch_w=TABLE_PATCH_W, int channels=TABLE_IMAGE_CHANNELS>
real patch_dist_exact(const Array<real> &a, const real *avec, const Array<real> &b, int bx, int by) {
#if TABLE_SSE
    __m128 ans(_mm_setzero_ps());
#else
    real ans = 0;
#endif
    const real *bvec = b.data + by * b.stride[0] + bx * b.stride[1];
    for (int dy = 0; dy < patch_w; dy++) {
        const real *arow = avec + dy * a.stride[0];
        const real *brow = bvec + dy * b.stride[0];
#if TABLE_SSE
        /* Compute 3*8 = 24 distances in 6 chunks of 4. TODO: This should be templated on the patch width */
#define PATCH_DIST_COMPUTE_DELTA(n) \
        __m128 a##n = _mm_loadu_ps(arow+n*4); \
        __m128 b##n = _mm_loadu_ps(brow+n*4); \
        __m128 delta##n = _mm_sub_ps(a##n, b##n); \
        delta##n = _mm_mul_ps(delta##n, delta##n);

        PATCH_DIST_COMPUTE_DELTA(0)
        PATCH_DIST_COMPUTE_DELTA(1)
        PATCH_DIST_COMPUTE_DELTA(2)
        PATCH_DIST_COMPUTE_DELTA(3)
        PATCH_DIST_COMPUTE_DELTA(4)
        PATCH_DIST_COMPUTE_DELTA(5)

        ans = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(ans, delta0), delta1), delta2), delta3), delta4), delta5);
#else
        for (int i = 0; i < TABLE_PATCH_W*channels; i++) {      // TODO The calling routine should be templated on the patch width
            real delta = arow[i]-brow[i];
            ans += delta*delta;
        }
#endif
//            if (early_terminate) {
//                if (ans > dmax) { return ans; }
//            }
    }
#if TABLE_SSE
    float *ans_f = (float *) &ans;
    return ans_f[0] + ans_f[1] + ans_f[2] + ans_f[3];
#else
    return ans;
#endif
}


template<class real, int TABLE_NCHANNELS>
real patch_dist_approx(const Array<real> &b_wh0, const real *avec, const real *bvec0, int xsrc, int ysrc) {
    const real *bvec = bvec0 + b_wh0.stride[0] * ysrc + b_wh0.stride[1] * xsrc;

#if !TABLE_SSE
    real dsrc = 0;
    for (int i = 0; i < TABLE_NCHANNELS; i++) {
#else
    /* Assume table dimensions is >= 4 and compute first 4 dimension distance in parallel */
    __m128 a0 = _mm_loadu_ps(avec);
    __m128 b0 = _mm_loadu_ps(bvec);
    __m128 delta = _mm_sub_ps(a0, b0);
    delta = _mm_mul_ps(delta, delta);

#define PATCH_DIST_APPROX_ADD_DELTA_BLOCK(n) \
        __m128 a##n = _mm_loadu_ps(avec+4*n); \
        __m128 b##n = _mm_loadu_ps(bvec+4*n); \
        __m128 delta##n = _mm_sub_ps(a##n, b##n); \
        delta##n = _mm_mul_ps(delta##n, delta##n); \
        delta = _mm_add_ps(delta, delta##n);
        
#define PATCH_DIST_APPROX_SUM_AND_RETURN() \
        float *deltaf = (float *) &delta; \
        return deltaf[0] + deltaf[1] + deltaf[2] + deltaf[3];

    if (TABLE_NCHANNELS == 8) {
        /* If table dimensions is exactly 8 then compute next 4 dimensions of distance also in parallel */
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(1);

        PATCH_DIST_APPROX_SUM_AND_RETURN();
    } else if (TABLE_NCHANNELS == 12) {
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(1);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(2);
        
        PATCH_DIST_APPROX_SUM_AND_RETURN();
    } else if (TABLE_NCHANNELS == 16) {
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(1);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(2);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(3);
        
        PATCH_DIST_APPROX_SUM_AND_RETURN();
    } else if (TABLE_NCHANNELS == 20) {
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(1);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(2);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(3);
        PATCH_DIST_APPROX_ADD_DELTA_BLOCK(4);

        PATCH_DIST_APPROX_SUM_AND_RETURN();
    } else {
        /* Otherwise, sum distance over remaining dimensions */
        float *deltaf = (float *) &delta;
        real dsrc = deltaf[0] + deltaf[1] + deltaf[2] + deltaf[3];
        for (int i = 4; i < TABLE_NCHANNELS; i++) {
#endif
            real delta = avec[i] - bvec[i];
            dsrc += delta*delta;
        }
        return dsrc;
#if TABLE_SSE
    }
#endif
}

#if TABLE_PATCHMATCH_EXACT
#define PATCHMATCH_PATCH_DIST(xsrc, ysrc) patch_dist_exact<real, TABLE_NCHANNELS>(a, avec, b, xsrc, ysrc)
#else
#define PATCHMATCH_PATCH_DIST(xsrc, ysrc) patch_dist_approx<real, TABLE_NCHANNELS>(wh_b, avec, bvec0, xsrc, ysrc)
#endif

#define REPLACE_IF_IMPROVED_MATCH() \
    { \
        real dcurrent = PATCHMATCH_PATCH_DIST(xsrc, ysrc); \
        if (dcurrent < dbest) { \
            dbest = dcurrent; \
            xbest = xsrc; \
            ybest = ysrc; \
        } \
    }
    
#if TABLE_DEBUG
#define PATCHMATCH_CHECK_KNN(check_next) \
    for (int check_k = 0; check_k < knn; check_k++) { \
        ASSERT(in_bounds(ann_L[check_k].x, bew), "expected existing bx in bounds"); \
        ASSERT(in_bounds(ann_L[check_k].y, beh), "expected existing by in bounds"); \
    } \
    if (x+1 < aew && (check_next) && (!allowed_patches_a || (*allowed_patches_a)(y, x+1) == p->allowed_index)) { \
        NNFPointer *ann_Lp = (NNFPointer *) &ann(y, x+1, 0, 0); \
        for (int check_k = 0; check_k < knn; check_k++) { \
            ASSERT(in_bounds(ann_Lp[check_k].x, bew), "expected existing bx at x+1 position in bounds"); \
            ASSERT(in_bounds(ann_Lp[check_k].y, beh), "expected existing by at x+1 position in bounds"); \
        } \
    }
#else
#define PATCHMATCH_CHECK_KNN(check_next)
#endif

#define REPLACE_IF_IMPROVED_MATCH_KNN() \
    { \
        if (min_dist2 > 0) { \
            int dx = xsrc - x; \
            int dy = ysrc - y; \
            int d = dx*dx+dy*dy; \
            if (d < min_dist2) { continue; } \
        } \
        \
        int cur_position = XY_TO_INT(xsrc, ysrc); \
        if (!prev_positions.count(cur_position)) { \
            real dcurrent = PATCHMATCH_PATCH_DIST(xsrc, ysrc); \
            if (dcurrent < ann_L[0].dist) { \
                int ann0_x = ann_L[0].x, ann0_y = ann_L[0].y; \
                prev_positions.erase(XY_TO_INT(ann0_x, ann0_y)); \
                prev_positions.insert(cur_position); \
                pop_heap(&ann_L[0], &ann_L[knn]); \
                ann_L[knn-1].x = xsrc; \
                ann_L[knn-1].y = ysrc; \
                ann_L[knn-1].dist = dcurrent; \
                push_heap(&ann_L[0],&ann_L[knn]); \
            } \
        } \
    }

#pragma pack(1)
class NNFPointer { public:
    double x, y, dist;
    
    INLINE bool operator < (const NNFPointer &b) const {
        return dist < b.dist;
    }
};

#define RS_SAMPLE_AROUND(xbest, ybest) \
    int xmin = xbest - rs_mag; \
    int xmax = xbest + rs_mag + 1; \
    int ymin = ybest - rs_mag; \
    int ymax = ybest + rs_mag + 1; \
    if (xmin < 0) { xmin = 0; } \
    if (ymin < 0) { ymin = 0; } \
    if (xmax > bew) { xmax = bew; } \
    if (ymax > beh) { ymax = beh; } \
    \
    int xsrc = xmin + rand() % (xmax-xmin); \
    int ysrc = ymin + rand() % (ymax-ymin);

template<class real, int is_descriptor, int TABLE_NCHANNELS, int use_knn=0, class itype=TABLE_DEFAULT_ITYPE>
void patchmatch(PatchTableParams *p, const Array<real> &a, const Array<real> &b, const Array<real> &wh_a, Array<real> &wh_b, Array<double> &ann, const Array<itype> *allowed_patches_a=NULL, const Array<itype> *allowed_patches_b=NULL, int astep=1) {
#if PM_VERBOSE
    printf("patchmatch\n"); fflush(stdout);
#endif
    double T0 = wall_time();
    int aew = gck_ew(wh_a), aeh = gck_eh(wh_a);
    int bew = gck_ew(wh_b), beh = gck_eh(wh_b);
    int knn = p->pm_knn;
    if (!use_knn) {
        ann.resize(aeh, aew, 3);
    } else {
        ann.resize(aeh, aew, knn, 3);
    }
    int min_dist2 = p->pm_min_dist*p->pm_min_dist;
    int prev_capacity = knn*8;
    int enrich_capacity = knn*knn*8;
    
    int aew_lo = (aew+astep-1)/astep;    // (aew_lo-1)*astep < aew
    int aeh_lo = (aeh+astep-1)/astep;

    if (allowed_patches_a) {
        if (allowed_patches_a->width() < aew || allowed_patches_a->height() < aeh) { fprintf(stderr, "allowed_patches_a size too small: %dx%d vs %dx%d\n", allowed_patches_a->width(), allowed_patches_a->height(), aew, aeh); ASSERT2(false, "allowed_patches_a incorrect size"); }
    }
    if (allowed_patches_b) {
        if (allowed_patches_b->width() < bew || allowed_patches_b->height() < beh) { fprintf(stderr, "allowed_patches_b size too small: %dx%d vs %dx%d\n", allowed_patches_b->width(), allowed_patches_b->height(), bew, beh); ASSERT2(false, "allowed_patches_b incorrect size"); }
    }

    int rs_max = MAX(bew, beh);

    bool do_enrich = (&wh_a == &wh_b) && p->pm_enrich;
#if PM_VERBOSE
    printf("patchmatch do_enrich=%d\n", int(do_enrich)); fflush(stdout);

    printf("patchmatch init\n"); fflush(stdout);
#endif
    PATCH_DIST_CALC_BVEC0(wh_b);
#if TABLE_OPENMP
    #pragma omp parallel
#endif
    {
        unordered_set<int> prev_positions(prev_capacity);
#if TABLE_OPENMP
        #pragma omp for
#endif
        for (int y = 0; y < aeh; y += astep) {
            int y_lo = y/astep;
            for (int x = 0; x < aew; x += astep) {
                int x_lo = x/astep;
                
                if (allowed_patches_a && ((*allowed_patches_a)(y, x) != p->allowed_index)) { continue; }
                
                PATCH_DIST_CALC_AVEC(wh_a);
                
                if (!use_knn) {
                    int bx = rand()%bew;
                    int by = rand()%beh;
                    ann(y_lo, x_lo, NNF_X) = bx;
                    ann(y_lo, x_lo, NNF_Y) = by;
                    ann(y_lo, x_lo, NNF_DIST) = PATCHMATCH_PATCH_DIST(bx, by);
                } else {
                    prev_positions.clear();
                    for (int k = 0; k < knn; k++) {
                        while (1) {
                            int bx = rand()%bew;
                            int by = rand()%beh;
                            if (allowed_patches_b && ((*allowed_patches_b)(by, bx) != p->allowed_index)) { continue; }
                            bool inserted = prev_positions.insert(XY_TO_INT(bx, by)).second;
                            if (!inserted) {
                                continue;
                            }
                            if (min_dist2 > 0) {
                                int dx = bx - x;
                                int dy = by - y;
                                int d = dx*dx+dy*dy;
                                if (d < min_dist2) { continue; }
                            }
                            ann(y_lo, x_lo, k, NNF_X) = bx;
                            ann(y_lo, x_lo, k, NNF_Y) = by;
                            ann(y_lo, x_lo, k, NNF_DIST) = PATCHMATCH_PATCH_DIST(bx, by);
                            break;
                        }
                    }
                    NNFPointer *ann_L = (NNFPointer *) &ann(y, x, 0, 0);
                    make_heap(ann_L, ann_L + knn);
                }
            }
        }
    }
    
    for (int iter = 0; iter < p->pm_iters; iter++) {
#if PM_VERBOSE
        printf("patchmatch iter %d/%d\n", iter, p->pm_iters); fflush(stdout);
#endif
        bool mirror = (iter % 2);
        int delta = mirror ? -1: 1;
        int bdelta = delta * astep;
#if TABLE_OPENMP
        for (int subtile = 0; subtile < 2; subtile++) {
#if PM_VERBOSE
            printf("patchmatch iter subtile %d\n", subtile); fflush(stdout);
#endif
            #pragma omp parallel
            {
#endif
                unordered_set<int> prev_positions(prev_capacity);
                unordered_set<int> enrich_positions(enrich_capacity);
#if TABLE_OPENMP
                int ithread = TABLE_OPENMP ? omp_get_thread_num(): 0;
                int nthreads = TABLE_OPENMP ? omp_get_num_threads(): 1;
                int itile = (ithread*2 + subtile);
                int ntiles = nthreads*2;
                int ytile_min = aeh*itile / ntiles;
                int ytile_max = aeh*(itile+1) / ntiles;
#else
                int ytile_min = 0, ytile_max = aeh;
#endif
                for (int y0 = ytile_min; y0 < ytile_max; y0++) {
                    int y = mirror ? (aeh - 1 - y0): y0;
                    if (y % astep != 0) { continue; }
                    int y_lo = y/astep;
                    
                    for (int x0 = 0; x0 < aew; x0++) {
                        int x = mirror ? (aew - 1 - x0): x0;
                        if (x % astep != 0) { continue; }
                        int x_lo = x/astep;
                        
                        if (allowed_patches_a && ((*allowed_patches_a)(y, x) != p->allowed_index)) { continue; }
                        
                        PATCH_DIST_CALC_AVEC(wh_a);
                        
                        if (!use_knn) {
                            int xbest = ann(y_lo, x_lo, NNF_X);
                            int ybest = ann(y_lo, x_lo, NNF_Y);
                            double dbest = ann(y_lo, x_lo, NNF_DIST);

                            /* Propagate x */
                            if (in_bounds(x_lo-delta, aew_lo) && (!allowed_patches_a || (*allowed_patches_a)(y, x-delta) == p->allowed_index)) {
                                int xsrc = ann(y_lo, x_lo-delta, NNF_X)+bdelta;
                                int ysrc = ann(y_lo, x_lo-delta, NNF_Y);
                                if (in_bounds(xsrc, bew) && (xsrc != xbest || ysrc != ybest) && (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index)) {
                                    REPLACE_IF_IMPROVED_MATCH();
                                }
                            }

                            /* Propagate y */
                            if (in_bounds(y_lo-delta, aeh_lo) && (!allowed_patches_a || (*allowed_patches_a)(y-delta, x) == p->allowed_index)) {
                                int xsrc = ann(y_lo-delta, x_lo, NNF_X);
                                int ysrc = ann(y_lo-delta, x_lo, NNF_Y)+bdelta;
                                if (in_bounds(ysrc, beh) && (xsrc != xbest || ysrc != ybest) && (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index)) {
                                    REPLACE_IF_IMPROVED_MATCH();
                                }
                            }
                            
                            /* Random search */
                            for (int rs_mag = rs_max; rs_mag >= 1; rs_mag /= 2) {
                                RS_SAMPLE_AROUND(xbest, ybest);
                                
                                if ((xsrc != xbest || ysrc != ybest) && (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index)) {
                                    REPLACE_IF_IMPROVED_MATCH();
                                }
                            }
                            
                            ann(y_lo, x_lo, NNF_X) = xbest;
                            ann(y_lo, x_lo, NNF_Y) = ybest;
                            ann(y_lo, x_lo, NNF_DIST) = dbest;
                        } else {
                            prev_positions.clear();
                            
                            NNFPointer *ann_L = (NNFPointer *) &ann(y_lo, x_lo, 0, 0);
                            for (int k = 0; k < knn; k++) {
                                int bx = ann_L[k].x;
                                int by = ann_L[k].y;
                                ASSERT(in_bounds(bx, bew), "expected existing bx in bounds");
                                ASSERT(in_bounds(by, beh), "expected existing by in bounds");
                                prev_positions.insert(XY_TO_INT(bx, by));
                            }
                            
                            /* Propagate x */
                            if (in_bounds(x_lo-delta, aew_lo) && (!allowed_patches_a || (*allowed_patches_a)(y, x-delta) == p->allowed_index)) {
                                for (int k = 0; k < knn; k++) {
                                    int xsrc = ann(y_lo, x_lo-delta, k, NNF_X)+bdelta;
                                    int ysrc = ann(y_lo, x_lo-delta, k, NNF_Y);

                                    ASSERT(in_bounds(ysrc, beh), "expected ysrc in bounds in prop x");
                                    if (in_bounds(xsrc, bew) && (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index)) {
                                        REPLACE_IF_IMPROVED_MATCH_KNN();
                                    }
                                }
                            }

                            /* Propagate y */
                            if (in_bounds(y_lo-delta, aeh_lo) && (!allowed_patches_a || (*allowed_patches_a)(y-delta, x) == p->allowed_index)) {
                                for (int k = 0; k < knn; k++) {
                                    int xsrc = ann(y_lo-delta, x_lo, k, NNF_X);
                                    int ysrc = ann(y_lo-delta, x_lo, k, NNF_Y)+bdelta;
                                    ASSERT(in_bounds(xsrc, bew), "expected xsrc in bounds in prop y");
                                    if (in_bounds(ysrc, beh) && (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index)) {
                                        REPLACE_IF_IMPROVED_MATCH_KNN();
                                    }
                                }
                            }

                            /* Random search */
                            for (int rs_mag = rs_max; rs_mag >= 1; rs_mag /= 2) {
                                for (int k = 0; k < knn; k++) {
                                    int xcurrent = ann_L[k].x;
                                    int ycurrent = ann_L[k].y;
                                    ASSERT(in_bounds(xcurrent, bew), "expected xcurrent in bounds in rs");
                                    ASSERT(in_bounds(ycurrent, beh), "expected ycurrent in bounds in rs");
                                    RS_SAMPLE_AROUND(xcurrent, ycurrent);
                                    ASSERT(in_bounds(xsrc, bew), "expected xsrc in bounds in rs");
                                    ASSERT(in_bounds(ysrc, beh), "expected ysrc in bounds in rs");
                                    if (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index) {
                                        REPLACE_IF_IMPROVED_MATCH_KNN();
                                    }
                                }
                            }

                            /* Forward enrichment */
                            if (do_enrich) {
                                enrich_positions.clear();
                                for (int k1 = 0; k1 < knn; k1++) {
                                    int x1 = ann_L[k1].x;
                                    int y1 = ann_L[k1].y;
                                    NNFPointer *ann_L2 = (NNFPointer *) &ann(y1/astep, x1/astep, 0, 0);
                                    for (int k2 = 0; k2 < knn; k2++) {
                                        int xsrc = ann_L2[k2].x;
                                        int ysrc = ann_L2[k2].y;
                                        ASSERT(in_bounds(xsrc, bew), "expected xcurrent in bounds in forward enrichment");
                                        ASSERT(in_bounds(ysrc, beh), "expected ycurrent in bounds in forward enrichment");
                                        if (!allowed_patches_b || (*allowed_patches_b)(ysrc, xsrc) == p->allowed_index) {
                                            bool inserted = enrich_positions.insert(XY_TO_INT(xsrc, ysrc)).second;
                                            if (inserted) {
                                                REPLACE_IF_IMPROVED_MATCH_KNN();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
#if TABLE_OPENMP
            }
        }
#endif
    }
    if (p->verbose) {
        printf("  patchmatch time: %f\n", wall_time()-T0); fflush(stdout);
    }
}

#define EMPTY_PATCH_INDEX -1
#define PATCH_INDEX(x, y) ((y)*(b0.width()-p->patch_w+1)+(x))
#define PATCH_INDEX_TO_X(i) ((i)%(b0.width()-p->patch_w+1))
#define PATCH_INDEX_TO_Y(i) ((i)/(b0.width()-p->patch_w+1))
    
#if TABLE_ENABLE_ANN
template<class real, class itype, int ndims>
shared_ptr<ANNkd_tree> build_ann_index_func(PatchTableParams *p, const Array<real> &wh0, const Array<itype> *allowed_patches, int step, vector<int> &ann_index_to_pos, ANNpointArray &ann_data_points) {
    int grid = step;
    int ntrain = DIV_ROUND_UP(gck_xmax(wh0)-gck_xmin(wh0), grid)*DIV_ROUND_UP(gck_ymax(wh0)-gck_ymin(wh0), grid);
    int ann_ntrain = ntrain;
    if (p->verbose) {
        printf("ann allocating data, %dx%d\n", ntrain, ndims);
    }
    ann_data_points = annAllocPts(ntrain, ndims);
    ann_index_to_pos.resize(ntrain);
    
    int train_index = 0;
    for (int y = gck_ymin(wh0); y < gck_ymax(wh0); y += grid) {
        for (int x = gck_xmin(wh0); x < gck_xmax(wh0); x += grid) {
            if (allowed_patches && (*allowed_patches)(y-gck_ymin(wh0), x-gck_xmin(wh0)) != p->allowed_index) { continue; }
            ANNpoint current_point = ann_data_points[train_index];
            for (int i = 0; i < ndims; i++) {
                current_point[i] = wh0(y, x, i);
            }
            ann_index_to_pos[train_index] = XY_TO_INT(x-gck_xmin(wh0), y-gck_ymin(wh0));
            train_index++;
        }
    }
    if (train_index > ntrain) {
        fprintf(stderr, "train_index (%d) > ntrain (%d), grid=%d\n", train_index, ntrain, grid); ASSERT2(false, "train_index != ntrain");
    }
    if (p->verbose) {
        printf("ann: creating tree\n");
    }
    
    shared_ptr<ANNkd_tree> ann_index(make_shared<ANNkd_tree>(ann_data_points, ntrain, ndims));
    
    if (p->verbose) {
        printf("ann: done building tree\n");
    }
    return ann_index;
}
    
#endif
    
#if TABLE_ENABLE_FLANN
typedef flann::Index<flann::L2<float> > FlannIndexType;
    
template<class real, class itype, int ndims>
shared_ptr<FlannIndexType> build_flann_index_func(PatchTableParams *p, const Array<real> &wh0, const Array<itype> *allowed_patches, int step, vector<int> &ann_index_to_pos) {
    double T0_build_index = wall_time();
    int ntrain = (gck_xmax(wh0)-gck_xmin(wh0))*(gck_ymax(wh0)-gck_ymin(wh0));
    
    if (p->verbose) {
        printf("build_flann_index: allocating train_data, %dx%d\n", ntrain, ndims);
    }
    float *train_data_ptr = new float[ntrain*ndims];
    ann_index_to_pos.resize(ntrain);
    
    int train_index = 0;
    TABLE_FOR_ALLOWED_PATCHES_STEP(step);
    
        for (int i = 0; i < ndims; i++) {
            //train_data[train_index][i] = wh0(y, x, i);
            train_data_ptr[train_index*ndims+i] = wh0(y, x, i);
        }
        ann_index_to_pos[train_index] = XY_TO_INT(x-gck_xmin(wh0), y-gck_ymin(wh0));
        
        train_index++;
    
    TABLE_END_FOR_ALLOWED_PATCHES();
    
    ntrain = train_index;
    flann::Matrix<float> train_data(train_data_ptr, ntrain, ndims);
    
    if (p->verbose) {
        printf("kdtree method: creating tree (ntrain=%d)\n", ntrain);
    }
    shared_ptr<FlannIndexType> flann_index(make_shared<FlannIndexType>(train_data, flann::KDTreeSingleIndexParams())); //, cv::cvflann::FLANN_DIST_L2);
    //            flann_index = make_shared<FlannIndexType>(train_data, flann::KDTreeIndexParams(p->flann_trees)); //, cv::cvflann::FLANN_DIST_L2);
    flann_index->buildIndex();
    
    delete[] train_data_ptr;
    if (p->verbose) {
        printf("kdtree method: done building tree (%f secs)\n", wall_time()-T0_build_index);
    }
    
    return flann_index;
}

#define TABLE_NCHANNELS ndims               /* Dimensions of descriptor (passed in manually or found by gck.h) */
#define GRID_NCHANNELS p->grid_ndims        /* First dimensions of descriptor are used for grid (GRID_NCHANNELS <= TABLE_NCHANNELS) */

template<class real, int is_descriptor, int ndims, int use_knn=0, class itype=TABLE_DEFAULT_ITYPE>
class TreeCANN { public:
    PatchTableParams *p;
    const Array<real> &b, wh_b;
    const Array<itype> *allowed_patches_b;
    int knn, matches_per_knn, spatial;
    vector<int> ann_index_to_pos;
#if TABLE_TREECANN_ANN
    shared_ptr<ANNkd_tree> index;
    ANNpointArray ann_data_points;
#else
    shared_ptr<FlannIndexType> index;
#endif
    
    TreeCANN(PatchTableParams *p_, const Array<real> &b_, const Array<real> &wh_b_, const Array<itype> *allowed_patches_b_=NULL, int knn_=1, int matches_per_knn_=4, int spatial_=1)
    :p(p_), b(b_), wh_b(wh_b_), allowed_patches_b(allowed_patches_b_), knn(knn_), matches_per_knn(matches_per_knn_), spatial(spatial_) {
        /* Build index */
#if TABLE_TREECANN_ANN
        index = build_ann_index_func<real, itype, ndims>(p, wh_b, allowed_patches_b, p->treecann_bgrid, ann_index_to_pos, ann_data_points);
#else
        index = build_flann_index_func<real, itype, ndims>(p, wh_b, allowed_patches_b, p->treecann_bgrid, ann_index_to_pos);
#endif
    }
    
    ~TreeCANN() {
#if TABLE_TREECANN_ANN
        if (ann_data_points) {
            annDeallocPts(ann_data_points);
            ann_data_points = NULL;
        }
#endif
    }

    void lookup(const Array<real> &a, const Array<real> &wh_a, Array<double> &ann, const Array<itype> *allowed_patches_a=NULL) {
        int agrid = p->treecann_agrid;
        int min_dist2 = p->pm_min_dist * p->pm_min_dist;
        
        int aew = gck_ew(wh_a), aeh = gck_eh(wh_a);
        int bew = gck_ew(wh_b), beh = gck_eh(wh_b);
        if (!use_knn) {
            ann.resize(aeh, aew, 3);
        } else {
            ann.resize(aeh, aew, knn, 3);
        }
        
        for (int y = 0; y < aeh; y++) {
            for (int x = 0; x < aew; x++) {
                if (allowed_patches_a && ((*allowed_patches_a)(y, x)) != p->allowed_index) { continue; }
                
                if (!use_knn) {
                    ann(y, x, NNF_X) = 0;
                    ann(y, x, NNF_Y) = 0;
                    ann(y, x, NNF_DIST) = DIST_INFINITY;
                } else {
                    for (int k = 0; k < knn; k++) {
                        ann(y, x, k, NNF_X) = 0;
                        ann(y, x, k, NNF_Y) = 0;
                        ann(y, x, k, NNF_DIST) = DIST_INFINITY;
                    }
                }
            }
        }

        int kmatch = knn*matches_per_knn;
        
#if TABLE_TREECANN_ANN
        ANNcoord query_vector[ndims];
        vector<ANNidx> matched_indices(kmatch);
        vector<ANNdist> matched_dists(kmatch);
#else
        flann::SearchParams search_params(p->flann_checks);
        search_params.eps = p->treecann_eps;
        
        vector<int> grid_index(ndims);
        vector<real> query_vector(ndims);
        
        flann::Matrix<real> query_matrix(&query_vector[0], 1, query_vector.size());
        vector<vector<int> > matched_indices_mat;
        vector<vector<float> > matched_dists_mat;
        
        matched_indices_mat.resize(1);
        matched_dists_mat.resize(1);
        matched_indices_mat[0].resize(kmatch);
        matched_dists_mat[0].resize(kmatch);
        
        ASSERT(matched_indices_mat.size() == 1, "expected 1 length for matched_indices_mat");
        ASSERT(matched_indices_mat[0].size() == kmatch, "expected kmatch length for matched_indices_mat[0]");
        vector<int> &matched_indices(matched_indices_mat[0]);
        vector<float> &matched_dists(matched_dists_mat[0]);
#endif
        vector<NNFPointer> matched_patch(kmatch);

        const Array<real> &wh0(wh_b);
        const Array<real> &b0(b);
        PATCH_DIST_CALC_BVEC0(wh0);
        
        unordered_set<int> prev_positions(knn*8);
        
        for (int y = 0; y < aeh; y += agrid) {
            for (int x = 0; x < aew; x += agrid) {
                if (allowed_patches_a && ((*allowed_patches_a)(y, x)) != p->allowed_index) { continue; }

                int x_wh = x+gck_xmin(wh_a);
                int y_wh = y+gck_ymin(wh_a);
                for (int i = 0; i < ndims; i++) {
                    query_vector[i] = wh_a(y_wh, x_wh, i);
                }
#if TABLE_TREECANN_ANN
                index->annkSearch(query_vector, kmatch, &matched_indices[0], &matched_dists[0], p->ann_eps);
#else
                index->knnSearch(query_matrix, matched_indices_mat, matched_dists_mat, kmatch, search_params);
#endif

                for (int k = 0; k < kmatch; k++) {
                    if (!in_bounds(matched_indices[k], ann_index_to_pos.size())) {
                        fprintf(stderr, "matched_indices[%d] = %d, out of bounds %d\n", k, matched_indices[k], int(ann_index_to_pos.size())); ASSERT2(false, "out of bounds");
                    }
                    matched_indices[k] = ann_index_to_pos[matched_indices[k]];
                }

                if (!use_knn) {
                    int xbest = ann(y, x, NNF_X);
                    int ybest = ann(y, x, NNF_Y);
                    double dbest = ann(y, x, NNF_DIST);
                    
                    {
                        PATCH_DIST_CALC_AVEC_EXACT_OR_APPROX(wh_a, a, x, y, p->is_descriptor);
                        
                        for (int k = 0; k < kmatch; k++) {
                            int v = matched_indices[k];
                            int xsrc = INT_TO_X(v);
                            int ysrc = INT_TO_Y(v);
                            double dcurrent = TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, p->is_descriptor);
                            if (dcurrent < dbest) {
                                dbest = dcurrent;
                                xbest = xsrc;
                                ybest = ysrc;
                            }
                        }
                    }
                    
                    for (int prop_y = -agrid; prop_y <= agrid; prop_y++) {
                        for (int prop_x = -agrid; prop_x <= agrid; prop_x++) {
                            int a_y_p = y + prop_y;
                            int a_x_p = x + prop_x;
                            if (allowed_patches_a && ((*allowed_patches_a)(a_y_p, a_x_p)) != p->allowed_index) { continue; }
                            if (!in_bounds(a_x_p, aew) || !in_bounds(a_y_p, aeh)) { continue; }
                            
                            PATCH_DIST_CALC_AVEC_EXACT_OR_APPROX(wh_a, a, a_x_p, a_y_p, p->is_descriptor);
                            for (int dy = -spatial; dy <= spatial; dy++) {
                                for (int dx = -spatial; dx <= spatial; dx++) {
                                    double dcurrent = 0;
                                    int xsrc = xbest + prop_x + dx;
                                    int ysrc = ybest + prop_y + dy;
                                    if (xsrc < 0) { xsrc = 0; }
                                    else if (xsrc >= bew) { xsrc = bew-1; }
                                    if (ysrc < 0) { ysrc = 0; }
                                    else if (ysrc >= beh) { ysrc = beh-1; }
                                    
                                    if (dx == 0 && dy == 0 && prop_x == 0 && prop_y == 0) {
                                        dcurrent = dbest;
                                    } else {
                                        dcurrent = TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, p->is_descriptor);
                                    }
                                    
                                    if (dcurrent < ann(a_y_p, a_x_p, NNF_DIST)) {
                                        ann(a_y_p, a_x_p, NNF_DIST) = dcurrent;
                                        ann(a_y_p, a_x_p, NNF_X) = xsrc;
                                        ann(a_y_p, a_x_p, NNF_Y) = ysrc;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    PATCH_DIST_CALC_AVEC_EXACT_OR_APPROX(wh_a, a, x, y, p->is_descriptor);

                    for (int k = 0; k < kmatch; k++) {
                        int v = matched_indices[k];
                        int xsrc = INT_TO_X(v);
                        int ysrc = INT_TO_Y(v);
                        matched_patch[k].dist = TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, p->is_descriptor);
                        matched_patch[k].x = xsrc;
                        matched_patch[k].y = ysrc;
                    }
                    
                    std::nth_element(matched_patch.begin(), matched_patch.begin()+knn-1, matched_patch.end());
                    for (int prop_y = -agrid; prop_y <= agrid; prop_y++) {
                        for (int prop_x = -agrid; prop_x <= agrid; prop_x++) {
                            int a_y_p = y + prop_y;
                            int a_x_p = x + prop_x;
                            if (allowed_patches_a && ((*allowed_patches_a)(a_y_p, a_x_p)) != p->allowed_index) { continue; }
                            if (!in_bounds(a_x_p, aew) || !in_bounds(a_y_p, aeh)) { continue; }

                            NNFPointer *ann_L = (NNFPointer *) &ann(a_y_p, a_x_p, 0, 0);
                            prev_positions.clear();
                            for (int k = 0; k < knn; k++) {
                                prev_positions.insert(XY_TO_INT(int(ann_L[k].x), int(ann_L[k].y)));
                            }

                            PATCH_DIST_CALC_AVEC_EXACT_OR_APPROX(wh_a, a, a_x_p, a_y_p, p->is_descriptor);
                            for (int k = 0; k < knn; k++) {
                                for (int dy = -spatial; dy <= spatial; dy++) {
                                    for (int dx = -spatial; dx <= spatial; dx++) {
                                        double dcurrent = 0;
                                        int xsrc = int(matched_patch[k].x) + prop_x + dx;
                                        int ysrc = int(matched_patch[k].y) + prop_y + dy;
                                        if (xsrc < 0) { xsrc = 0; }
                                        else if (xsrc >= bew) { xsrc = bew-1; }
                                        if (ysrc < 0) { ysrc = 0; }
                                        else if (ysrc >= beh) { ysrc = beh-1; }
                                        
                                        if (min_dist2 > 0) {
                                            int dx = xsrc - a_x_p;
                                            int dy = ysrc - a_y_p;
                                            int d = dx*dx+dy*dy;
                                            if (d < min_dist2) { continue; }
                                        }
                                        
                                        int cur_position = XY_TO_INT(xsrc, ysrc);
                                        if (!prev_positions.count(cur_position)) {
                                            real dcurrent = TABLE_LOOKUP_PATCH_DIST_IS_DESCRIPTOR(xsrc, ysrc, p->is_descriptor);
                                            if (dcurrent < ann_L[0].dist) {
                                                int ann0_x = ann_L[0].x, ann0_y = ann_L[0].y;
                                                prev_positions.erase(XY_TO_INT(ann0_x, ann0_y));
                                                prev_positions.insert(cur_position); 
                                                pop_heap(&ann_L[0], &ann_L[knn]); 
                                                ann_L[knn-1].x = xsrc; 
                                                ann_L[knn-1].y = ysrc; 
                                                ann_L[knn-1].dist = dcurrent; 
                                                push_heap(&ann_L[0],&ann_L[knn]); 
                                            } 
                                        } 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};
#endif

/* ------------------------------------------------------------------------------------------------
   Partitions -- divide grid cells non-uniformly along each dimension
   ------------------------------------------------------------------------------------------------ */

template<class real>
class PatchComparator { public:
    const Array<real> &wh0;
    int channel;
    PatchComparator(const Array<real> &wh0_, int channel_) :wh0(wh0_) {
        channel = channel_;
    }
    INLINE bool operator ()(int a, int b) {
        int ax = INT_TO_X(a);
        int ay = INT_TO_Y(a);
        int bx = INT_TO_X(b);
        int by = INT_TO_Y(b);
        return wh0(ay, ax, channel) < wh0(by, bx, channel);
    }
};

template<class real>
class PatchPartition1D { public:
    vector<int> to_bin;             /* Remap int index in [0, PATCH_PARTITION_NMAP) to a bin index */
    vector<real> centers;           /* Center of each cluster */
    vector<real> lsizes;            /* Half size of each cluster in the negative (left) direction */
    vector<real> rsizes;            /* Half size of each cluster in the positive (right) direction */
    real min_val, max_val;
    real range;
    int nslices;
    PatchTableParams *p;
#if TABLE_ALLOW_COMPARE
    vector<double> to_bin_frac;
#endif

    void save(FILE *f) {
        ::save(f, "PatchPartition1D");
        ::save(f, to_bin);
        ::save(f, centers);
        ::save(f, lsizes);
        ::save(f, rsizes);
        ::save(f, min_val);
        ::save(f, max_val);
        ::save(f, range);
        ::save(f, nslices);
        ::save(f, "EndPatchPartition1D");
    }
    
    PatchPartition1D(PatchTableParams *p_, FILE *f) {
        p = p_;
        load_id(f, "PatchPartition1D");
        load(f, to_bin);
        load(f, centers);
        load(f, lsizes);
        load(f, rsizes);
        load(f, min_val);
        load(f, max_val);
        load(f, range);
        load(f, nslices);
        load_id(f, "EndPatchPartition1D");
    }
    
    PatchPartition1D(PatchTableParams *p_, const Array<real> &wh0, int channel, int nslices_, vector<int> patches, bool randomize, const vector<real> &min_valL, const vector<real> &max_valL) {
        p = p_;
        nslices = nslices_;
        
        ASSERT2(nslices >= 1, "expected nslices >= 1");
        sort(patches.begin(), patches.end(), PatchComparator<real>(wh0, channel));
        min_val = min_valL[channel]; //wh0(INT_TO_Y(patches[0]), INT_TO_X(patches[0]), channel);
        max_val = max_valL[channel]; //wh0(INT_TO_Y(patches[patches.size()-1]), INT_TO_X(patches[patches.size()-1]), channel);
        range = max_val - min_val;
        
        centers.resize(nslices);
        lsizes.resize(nslices);
        rsizes.resize(nslices);
        if (p->regular_grid) {                      /* TODO: Optimization: in regular_grid mode don't need to sort */
            double bin_size = (max_val-min_val)*1.0/nslices;
            for (int i = 0; i < nslices; i++) {
                double t = (i+0.5)/(nslices);
                centers[i] = min_val + (max_val-min_val)*t;
                lsizes[i] = bin_size/2;
                rsizes[i] = bin_size/2;
            }
        } else {
            for (int i = 0; i < nslices; i++) {
                int idx1 = int(patches.size()-1)*i/nslices;
                int idx2 = int(patches.size()-1)*(i+1)/nslices;
                real cell_max = wh0(INT_TO_Y(patches[idx2]), INT_TO_X(patches[idx2]), channel);
                real cell_min = wh0(INT_TO_Y(patches[idx1]), INT_TO_X(patches[idx1]), channel);
                centers[i] = (cell_min+cell_max)/2;
                if (p->verbose && i == 0) { printf("centers[%d] = %f, idx1=%d, idx2=%d, cell_min=%f, cell_max=%f, patches.size()=%d\n", i, centers[i], idx1, idx2, cell_min, cell_max, patches.size()); }
            }
            
            for (int i = 0; i < nslices; i++) {
                real prev_center = (i-1 >= 0) ? centers[i-1]: min_val;
                real next_center = (i+1 < (int) centers.size()) ? centers[i+1]: max_val;
                lsizes[i] = (centers[i]-prev_center);
                if (i-1 >= 0) { lsizes[i] /= 2; }
                rsizes[i] = (next_center-centers[i]);
                if (i+1 < (int) centers.size()) { rsizes[i] /= 2; }
                if (lsizes[i] < 0) { fprintf(stderr, "lsizes[%d] = %f, centers[i]=%f, prev_center=%f, min_val=%f\n", i, lsizes[i], centers[i], prev_center, min_val); ASSERT2(false, "expected lsizes >= 0"); }
                if (rsizes[i] < 0) { fprintf(stderr, "rsizes[%d] = %f\n", i, rsizes[i]); ASSERT2(false, "expected rsizes >= 0"); }
#if TABLE_EXTRA_VERBOSE
                printf("channel %d, i=%d, nslices=%d, patches=%d, centers=%f, lsizes=%f, rsizes=%f\n", channel, i, nslices, patches.size(), centers[i], lsizes[i], rsizes[i]);
#endif
            }
        }
        
        to_bin.resize(PATCH_PARTITION_NMAP);
#if TABLE_ALLOW_COMPARE
        to_bin_frac.resize(PATCH_PARTITION_NMAP);
#endif
        int icenter = 0;
        for (int i = 0; i < PATCH_PARTITION_NMAP; i++) {
            double coord = min_val + (i+0.5) * range / (PATCH_PARTITION_NMAP);
            while (coord > centers[icenter]) {       /* Find the smallest center to the right of us */
                if (icenter+1 < int(centers.size())) { icenter++; }
                else { break; }
            }
            int prev_center = icenter-1;
            if (prev_center < 0) { prev_center = 0; }
            int next_center = icenter;
            if (next_center == 0) { next_center++; }
            if (nslices == 1) { next_center = prev_center = 0; }
            double mid = (centers[prev_center] + centers[next_center]) / 2.0;
            
            //        mid01    mid12
            //          |        |
            //     c0       c1       c2
            // ----*--------*--------*----
            //
            //int bin_idx = i * nslices / PATCH_PARTITION_NMAP;

            int bin_idx;
            if (!randomize) {
                bin_idx = (coord < mid) ? prev_center: next_center;
            } else {
                double prob = (coord - centers[prev_center]) / (centers[next_center] - centers[prev_center]);
                bin_idx = (rand_f() < prob) ? next_center: prev_center;
            }
            ASSERT(in_bounds(bin_idx, nslices), "expected bin_idx in bounds");
            to_bin[i] = bin_idx;

#if TABLE_ALLOW_COMPARE
            double bin_idx_frac = 0.0;
            if (coord < mid) {
                double t_frac = (coord - centers[prev_center]) / (mid - centers[prev_center]);
                to_bin_frac[i] = prev_center + t_frac * 0.5;
            } else {
                double t_frac = (coord - mid) / (centers[next_center] - mid);
                to_bin_frac[i] = prev_center + 0.5 + t_frac * 0.5;
            }
#endif

#if TABLE_EXTRA_VERBOSE
            if (i == 588) {
                printf("PatchPartition channel %d, i=%d, coord=%f, prev_center=%d (%f), next_center=%d (%f), mid=%f, bin_idx=%d\n", channel, i, coord, prev_center, centers[prev_center], next_center, centers[next_center], mid, bin_idx);
            }
#endif
        }
        
#if TABLE_VERBOSE
        printf("PatchPartition1D(channel=%d)\n", channel);
        printf("  min: %f, max: %f\n", double(min_val), double(max_val));
        printf("  centers: ");
        for (int i = 0; i < nslices; i++) {
            printf("%f, ", double(centers[i]));
        }
        printf("\n");
        printf("  mapping: ");
        for (int i = 0; i < PATCH_PARTITION_NMAP; i += PATCH_PARTITION_NMAP/32) {
            printf("%d, ", to_bin[i]);
        }
        printf("\n");
        printf("\n");
#endif
    }
    
#if TABLE_ALLOW_COMPARE
    INLINE double get_bin_frac(real v) {
        int idx = (v-min_val)*(PATCH_PARTITION_NMAP-1)/range;
        if (idx < 0) { idx = 0; }
        else if (idx >= PATCH_PARTITION_NMAP) { idx = PATCH_PARTITION_NMAP-1; }
        ASSERT(in_bounds(idx, PATCH_PARTITION_NMAP), "expected get_bin() index in range");
        return to_bin_frac[idx];
    }
#endif

    INLINE int get_bin(real v) {
        int idx = (v-min_val)*(PATCH_PARTITION_NMAP-1)/range;
        if (idx < 0) { idx = 0; }
        else if (idx >= PATCH_PARTITION_NMAP) { idx = PATCH_PARTITION_NMAP-1; }
        ASSERT(in_bounds(idx, PATCH_PARTITION_NMAP), "expected get_bin() index in range");
        int ans = to_bin[idx];
        /*
        // This commented-out code performs some sanity checks on the grids
        int ans_brute = -1;
        double dbrute = 1e100;
        for (int i = 0; i < nslices; i++) {
            double dcurrent = fabs(v-centers[i]);
            if (dcurrent < dbrute) {
                dbrute = dcurrent;
                ans_brute = i;
            }
        }
        int ans_next = idx+1 < PATCH_PARTITION_NMAP ? to_bin[idx+1]: -1;
        int ans_prev = idx-1 >= 0 ? to_bin[idx-1]: -1;
        printf("get_bin(%f) => %d (idx: %d, next: %d, prev: %d), (%f), %d (%f)\n", double(v), ans, idx, ans_next, ans_prev, fabs(centers[ans]-v), ans_brute, fabs(centers[ans_brute]-v));
        if (ans != ans_brute) {
            if (!(ans_prev <= ans_brute && ans_next >= ans_brute)) {
                printf("get_bin ans != ans_brute\n"); exit(1);
            }
        }*/
        
        return ans;
    }
};

/* Finds bin partitions along all dimensions */
template<class real, class itype>
class PatchPartition { public:
    vector<PatchPartition1D<real> *> L;
    PatchTableParams *p;
    
    void save(FILE *f) {
        ::save(f, "PatchPartition");
        ::save(f, int(L.size()));
        for (int i = 0; i < (int) L.size(); i++) {
            L[i]->save(f);
        }
        ::save(f, "EndPatchPartition");
    }
    
    PatchPartition(PatchTableParams *p_, FILE *f) {
        p = p_;
        load_id(f, "PatchPartition");
        int size;
        load(f, size);
        L.resize(size);
        for (int i = 0; i < (int) L.size(); i++) {
            L[i] = new PatchPartition1D<real>(p_, f);
        }
        load_id(f, "EndPatchPartition");
    }

    PatchPartition(PatchTableParams *p_, const Array<real> &wh0, const vector<int> &nslices, bool randomize, const Array<itype> *allowed_patches, const vector<real> &min_valL, const vector<real> &max_valL) {
        p = p_;
        vector<int> patches;
        int h = (gck_ymax(wh0) - gck_ymin(wh0) + p->partition_step-1) / p->partition_step;
        int w = (gck_xmax(wh0) - gck_xmin(wh0) + p->partition_step-1) / p->partition_step;
        patches.resize(w*h);
        int count = 0;
        
        if (!allowed_patches) {
            for (int y = gck_ymin(wh0); y < gck_ymax(wh0); y += p->partition_step) {
                for (int x = gck_xmin(wh0); x < gck_xmax(wh0); x += p->partition_step) {
                    patches[count++] = XY_TO_INT(x, y);
                }
            }
        } else {
            for (int y = gck_ymin(wh0); y < gck_ymax(wh0); y += p->partition_step) {
                for (int x = gck_xmin(wh0); x < gck_xmax(wh0); x += p->partition_step) {
                    if ((*allowed_patches)(y-gck_ymin(wh0), x-gck_xmin(wh0)) == p->allowed_index) {
                        patches[count++] = XY_TO_INT(x, y);
                    }
                }
            }
        }
        if (!allowed_patches) {
            ASSERT2(count == w*h, "expected patch count to be w*h");
        } else {
            patches.resize(count);
        }
        
        for (int i = 0; i < (int) nslices.size(); i++) {
            L.push_back(new PatchPartition1D<real>(p, wh0, i, nslices[i], patches, randomize, min_valL, max_valL));
        }
    }
    
    ~PatchPartition() {
        for (int i = 0; i < (int) L.size(); i++) {
            delete L[i];
        }
    }
    
    void get_patch_index(const Array<real> &wh, int x, int y, vector<int> &patch_index) {
        ASSERT(patch_index.size() == L.size(), "expected patch_index of same size as L");
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            patch_index[i] = L[i]->get_bin(wh(y, x, i));
        }
    }
    
    template<int use_dmax>
    float patch_dist_to_grid(const Array<real> &wh, int x, int y, const vector<int> &grid_index, float dmax=0, bool verbose=false) {
        /* TODO: Could also add a penalty that sums over all WH dimensions (including those not added to the table), patch distance from their mean */
        ASSERT(grid_index.size() == GRID_NCHANNELS, "expected grid_index size to equal table channels");
        float ans = 0;
        for (int i = 0; i < (int) grid_index.size(); i++) {
            ASSERT(in_bounds(grid_index[i], L[i]->centers.size()), "expected grid_index to be in range");
            float delta = wh(y, x, i) - L[i]->centers[grid_index[i]];
#if TABLE_INFINITY_NORM
            float sz = delta > 0 ? L[i]->rsizes[grid_index[i]]: L[i]->lsizes[grid_index[i]];
            if (verbose) {
                printf("patch_dist_to_grid: wh=%f, delta[%d]=%f, sz=%f, |delta|/sz=%f, grid_index=%d, centers=%f\n", wh(y, x, i), i, delta, sz, fabs(delta)/sz, grid_index[i], L[i]->centers[grid_index[i]]);
            }
            delta = fabs(delta) / sz;
            ans = MAX(ans, delta);
#else
            ans += delta*delta;
#endif
            if (use_dmax && ans > dmax) { return ans; }
        }
        if (verbose) {
            printf("patch_dist_to_grid: return %f\n\n", ans);
        }
        return ans;
    }
};

/* Tracks for integers i = 0...n-1, a set of integers j = 0...n-1 that are adjacent to i. */
class AdjacencySet { public:
    bool unique;
    vector<unordered_set<int> *> sets;
    vector<vector<int> *> setsL;
    
    AdjacencySet(int n, bool unique_);
    ~AdjacencySet();
    void add(int i, int j);
    void compute_sets();

    INLINE vector<int> *get_set(int i) {
        return setsL[i];
    }
};

vector<int> scale_slices_to_limit(vector<int> nslices, double limit);

class PropEdge { public:
    int dim;
    int dir;
    int label;
    PropEdge(int dim_, int dir_, int nlabels) :dim(dim_), dir(dir_) {
        label = (2*(dim)+1+((dir)+1)/2);
        if (!in_bounds(label, nlabels)) {
            fprintf(stderr, "label out of bounds in PropEdge constructor: dim=%d, dir=%d, label=%d, nlabels=%d\n", dim, dir, label, nlabels);
            ASSERT2(in_bounds(label, nlabels), "label out of bounds in PropEdge constructor");
        }
    }
};

#pragma pack(1)
template<class itype>
class PropElement { public:
    itype grid_index;
    itype table_value;
    uint8_t label;
    PropElement() { }
    PropElement(itype grid_index_, itype table_value_, itype label_=0) :grid_index(grid_index_), table_value(table_value_), label(label_) { }
};

template<class itype>
class PropElementProductQuantize { public:
    itype grid_index;
    itype table_value;
    float path_length;
    PropElementProductQuantize() { }
    PropElementProductQuantize(itype grid_index_, itype table_value_, float path_length_) :grid_index(grid_index_), table_value(table_value_), path_length(path_length_) { }
    
    INLINE bool operator < (const PropElementProductQuantize &b) const {
        return path_length > b.path_length;         /* Use > comparison to create a min heap */
    }
};

/* ------------------------------------------------------------------------------------------------
   PatchTable class with compile-time dimension count
   ------------------------------------------------------------------------------------------------ */
    
template<class real, class in_type, int ndims, class itype=TABLE_DEFAULT_ITYPE>
class PatchTableFixedN {
  private:
#if TABLE_MULTI_TABLES
    vector<shared_ptr<PatchTableFixedN<real, in_type, ndims, itype> > > tables;     /* List of (ntables-1) other tables */
    vector<shared_ptr<PatchTableParams> > tables_params;
#endif
    
#if TABLE_PRODUCT_QUANTIZE
    shared_ptr<ProductQuantizer<real, itype> > product_quantizer;
#endif
    
#if TABLE_CLUSTER_KMEANS
    vector<shared_ptr<PatchTableFixedN<real, in_type, ndims, itype> > > sub_tables;  /* Sub-tables, one for each k-means cluster */
    Array<itype> sub_allowed_patches;
    int sub_upper_left;
#endif

    AdjacencySet *table_sets;                               /* Only used internally by constructor */
    vector<real> min_val, max_val, range;                   /* Only used internally by constructor */
    Array<in_type> lookup_buffer;                           /* Only used internally by lookup() method: stores YUV color image */
    Array<real> lookup_wh;                                  /* Only used internally by lookup() method: stores WH (GCK) coordinates */
    Array<real> lookup_dist_orig;                           /* Only used internally by lookup() method */
    shared_ptr<Array<itype> > allowed_patches;              /* Used internally by constructor and lookup method */

#if TABLE_KCOHERENCE_TRIANGLE
    typedef pair<int, float> kcoherence_type;
#define TABLE_KCOHERENCE_POS(v) ((v).first)
#define TABLE_KCOHERENCE_DIST(v) ((v).second)
#else
    typedef int kcoherence_type;
#define TABLE_KCOHERENCE_POS(v) (v)
#endif
    
    Array<kcoherence_type> kcoherence_set;
    shared_ptr<cv::PCA> pca;
#if TABLE_ENABLE_ANN
    shared_ptr<ANNkd_tree> ann_index;
    ANNpointArray ann_data_points;
#endif

#if TABLE_ENABLE_FLANN
    shared_ptr<FlannIndexType> flann_index;
    vector<int> kdtree_index_to_pos;
    shared_ptr<TreeCANN<real, 1, ndims, 0, itype> > treecann;
#endif

#if (TABLE_ENABLE_ANN || TABLE_ENABLE_FLANN)
    vector<int> ann_index_to_pos;
#endif
    
#if (TABLE_OPTIMIZE_DT && !TABLE_DT_REGULAR)
    Array<dist_t> table_dist;                               /* Only used internally by constructor */
#endif
    int table_add_count;
    
  public:
    Array<itype> table;

    PatchTableParams *p;

    Array<in_type> b0;                                      /* Database image. We find a NNF from a => b. */
    Array<real> wh0;                                        /* Walsh-Hadamard (Grey-Code Kernels) transform of database image */
    vector<int> nslices;
    PatchPartition<real, itype> *part;

#if TABLE_DEBUG
    void check_kcoherence() {
        if (p->kcoherence_step == 1) {
            printf("check_kcoherence\n");
            for (int y = 0; y < kcoherence_set.height(); y++) {
                for (int x = 0; x < kcoherence_set.width(); x++) {
                    if (allowed_patches && (*allowed_patches)(y, x) != p->allowed_index) { continue; }
                    for (int k = 0; k < kcoherence_set.channels(); k++) {
                        int v = TABLE_KCOHERENCE_POS(kcoherence_set(y, x, k));
                        int xp = INT_TO_X(v);
                        int yp = INT_TO_Y(v);
                        if (!in_bounds(yp, b0.height()-p->patch_w+1) || !in_bounds(xp, b0.width()-p->patch_w+1)) {
                            fprintf(stderr, "x,y,k=%d,%d,%d, xp,yp=%d,%d out of bounds, b0 width: %d, height: %d, patch_w=%d\n", x, y, k, xp, yp, b0.width(), b0.height()); ASSERT2(false, "out of bounds");
                        }
                    }
                }
            }
        }
    }
#endif

    void save(FILE *f) {
        ::save(f, "PatchTable");
        ::save(f, table);
        ::save(f, b0);
        ::save(f, wh0);
        ::save(f, nslices);
        part->save(f);
        ::save(f, "EndPatchTable");
    }
    
    void load(PatchTableParams *p_, FILE *f) {
        p = p_;
        ::load_id(f, "PatchTable");
        ::load(f, table);
        ::load(f, b0);
        ::load(f, wh0);
        ::load(f, nslices);
        part = new PatchPartition<real, itype>(p_, f);
        load_id(f, "EndPatchTable");
    }

    void save(string filename) {
        double T0 = wall_time();
        FILE *f = fopen(filename.c_str(), "wb");
        save(f);
        fclose(f);
        if (p->verbose) {
            printf("PatchTable saved to %s in %f secs\n", filename.c_str(), wall_time()-T0);
        }
    }
    
    void load(PatchTableParams *p_, string filename) {
        p = p_;
        double T0 = wall_time();
        FILE *f = fopen(filename.c_str(), "rb");
        if (!f) { fprintf(stderr, "PatchTable: could not load from %s\n", filename.c_str()); exit(1); }
        load(p_, f);
        fclose(f);
        if (p->verbose) {
            printf("PatchTable loaded from %s in %f secs\n", filename.c_str(), wall_time()-T0);
        }
    }
    
    PatchTableFixedN(PatchTableParams *p_, string filename) {
#if TABLE_ENABLE_ANN
        ann_data_points = NULL;
#endif
        load(p_, filename);
    }

    void add_patch(int x, int y, const vector<int> &grid_index, PatchTableParams *p, bool unique) {
#if (TABLE_DEBUG && 0)
        for (int i = 0; i < (int) grid_index.size(); i++) {
            if (grid_index[i] < 0) { fprintf(stderr, "grid_index[%d] = %d, out of bounds, below zero\n", i, grid_index[i]); exit(1); }
            if (grid_index[i] >= nslices[i]) { fprintf(stderr, "grid_index[%d] = %d, out of bounds, above slices\n", i, grid_index[i]); exit(1); }
        }
#endif
#if TABLE_EXTRA_VERBOSE
        if (!in_bounds(x-gck_xmin(wh0), b0.width()-p->patch_w+1) || !in_bounds(y-gck_ymin(wh0), b0.height()-p->patch_w+1)) { fprintf(stderr, "add_patch source location out of bounds: %d, %d, %dx%d, patch_w=%d\n", x-gck_xmin(wh0), y-gck_ymin(wh0), b0.width()-p->patch_w+1, b0.height()-p->patch_w+1, p->patch_w); ASSERT2(false, "out of bounds"); exit(1); }
#endif
        int v = table(grid_index);
        if (v == TABLE_UNUSED) {
            table(grid_index) = XY_TO_INT(x, y);
#if !TABLE_DT_REGULAR
#if TABLE_OPTIMIZE_DT
#if TABLE_OPTIMIZE_REGULAR
            table_dist(grid_index) = 0;
#else
            table_dist(grid_index) = part->template patch_dist_to_grid<0>(wh0, x, y, grid_index);
#endif
#endif
#endif
#if TABLE_DT_REGULAR
            table_add_count++;
#endif
        } else {
#if !TABLE_RANDOMIZE_ADD
            int x_orig = INT_TO_X(v);
            int y_orig = INT_TO_Y(v);
            if (x == x_orig && y == y_orig) { return; }
            
#if (!TABLE_OPTIMIZE_DT || TABLE_OPTIMIZE_REGULAR || TABLE_DT_REGULAR)
            float d_orig = part->template patch_dist_to_grid<0>(wh0, x_orig, y_orig, grid_index);
#else
            dist_t d_orig = table_dist(grid_index);
#endif
            float d_ours = part->template patch_dist_to_grid<1>(wh0, x, y, grid_index, d_orig);
            if (p->randomize_dt) {
                int i1 = PATCH_INDEX(x-gck_xmin(wh0), y-gck_ymin(wh0));
                int i2 = PATCH_INDEX(x_orig-gck_xmin(wh0), y_orig-gck_ymin(wh0));
                if (!unique) {
                    table_sets->add(i1, i2);        // TODO: In the new unique mode this shouldn't be needed
                }
                table_sets->add(i2, i1);
            }
            if (d_ours < d_orig && !unique) {
                table(grid_index) = XY_TO_INT(x, y);
#if (TABLE_OPTIMIZE_DT && !TABLE_OPTIMIZE_REGULAR && !TABLE_DT_REGULAR)
                table_dist(grid_index) = d_ours;
#endif
            }
#else       /* Case TABLE_RANDOMIZE_ADD */
            if (rand()%2 == 0) {
                dist_t d_ours = part->template patch_dist_to_grid<0>(wh0, x, y, grid_index);
                table(grid_index) = XY_TO_INT(x, y);
#if TABLE_OPTIMIZE_DT
                table_dist(grid_index) = d_ours;
#endif
            }
#endif
        }
    }
    
    void reduce_dim(const Array<in_type> &a, Array<real> &wh) {
        bool do_filter = p->filter_dims.size();
        
        int ndims0 = p->ndims;
        if (do_filter) {
            if (p->filter_dims.size() != p->ndims) {
                fprintf(stderr, "expected filter_dims size to match p->ndims\n");
                ASSERT2(false, "filter_dims size != p->ndims");
            }
            int max_filter_dim = max<int>(p->filter_dims);
            p->ndims = max_filter_dim+1;
            if (p->dim_algo == TABLE_DIM_ALGO_WH) {
                if (p->ndims > 30) { p->ndims = 40; }       // TODO: This should be generated from supported dimensions automatically
                else if (p->ndims > 20) { p->ndims = 30; }
                else if (p->ndims > 12) { p->ndims = 20; }
            }
        }
        
        Array<real> wh_temp;
        
        if (p->dim_algo == TABLE_DIM_ALGO_WH) {
            gck<in_type, real>(a, do_filter ? wh_temp: wh, p->ndims-p->nchroma*2, p->nchroma, p->patch_w);
        } else if (p->dim_algo == TABLE_DIM_ALGO_PCA) {
            if (!pca) {
                pca = get_patch_pca<in_type>(p, a);
            }
            apply_patch_pca<in_type>(p, a, pca, do_filter ? wh_temp: wh);
        }
        
        if (do_filter) {
            wh.resize(wh_temp.height(), wh_temp.width(), ndims0);
            for (int y = gck_ymin(wh); y < gck_ymax(wh); y++) {
                for (int x = gck_xmin(wh); x < gck_xmax(wh); x++) {
                    for (int k = 0; k < ndims0; k++) {
                        wh(y, x, k) = wh_temp(y, x, p->filter_dims[k]);
                    }
                }
            }
            
            p->ndims = ndims0;
        }
    }
    
    void count_table(string info="table count after populate") {
#if (TABLE_DEBUG && 0)
        if (p->verbose) {
            printf("begin count_table\n");
            vector<int> patch_used((b0.width()-p->patch_w+1)*(b0.height()-p->patch_w+1), 0);
            for (int j = 0; j < table.nelems; j++) {
                int v0 = table.data[j];
                if (v0 != TABLE_UNUSED) {
                    int xsrc = INT_TO_X(v0)-gck_xmin(wh0);
                    int ysrc = INT_TO_Y(v0)-gck_ymin(wh0);
                    ASSERT2((unsigned) xsrc < (unsigned) (gck_xmax(wh0) - gck_xmin(wh0)) &&
                            (unsigned) ysrc < (unsigned) (gck_ymax(wh0) - gck_ymin(wh0)), "expected table patch to be in bounds");
                    int isrc = PATCH_INDEX(xsrc, ysrc);
                    ASSERT2((unsigned) isrc < (unsigned) patch_used.size(), "expected isrc in range");
                    patch_used[isrc] = 1;
                }
            }
            int table_count = 0;
            for (int i = 0; i < (int) patch_used.size(); i++) {
                if (patch_used[i]) { table_count++; }
            }
            int total_count = (b0.height()-p->patch_w+1)*(b0.width()-p->patch_w+1);
            printf("  table cells: %d\n", table.nelems);
            printf("  %s: %d (%.1f%% of %d possible patches)\n", info.c_str(), table_count, table_count*100.0/total_count, total_count);
        }
#endif
    }
    
    void set_min_max(const Array<real> &choose_wh0) {
        min_val.resize(GRID_NCHANNELS);
        max_val.resize(GRID_NCHANNELS);
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            min_val[i] = 1e20;
            max_val[i] = -1e20; //max_val[i] = wh0(0, 0, i);
        }
        
        TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER()
            for (int i = 0; i < GRID_NCHANNELS; i++) {
                real current_val = choose_wh0(y, x, i);
                min_val[i] = MIN(min_val[i], current_val);
                max_val[i] = MAX(max_val[i], current_val);
            }
        TABLE_END_FOR_ALLOWED_PATCHES()
        
        range.resize(GRID_NCHANNELS);
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            range[i] = max_val[i]-min_val[i];
        }
    }
    
    void set_default_nslices(const Array<real> &choose_wh0) {
        nslices.resize(GRID_NCHANNELS);

        vector<real> sizeL(range);
#if TABLE_SIZE_STDDEV
        vector<real> sum_x(GRID_NCHANNELS), sum_x2(GRID_NCHANNELS);
        int patch_count = 0;
        TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER()
            for (int i = 0; i < GRID_NCHANNELS; i++) {
                real current_val = choose_wh0(y, x, i);
                sum_x[i] += current_val;
                sum_x2[i] += current_val*current_val;
            }
            patch_count++;
        TABLE_END_FOR_ALLOWED_PATCHES()

        for (int i = 0; i < GRID_NCHANNELS; i++) {
            double E = sum_x[i] / patch_count;
            double V = sum_x2[i]  / patch_count - E*E;
            if (V < 0) { fprintf(stderr, "Warning: clamping variance less than zero: %f", V); V = 0; }
            double sigma = sqrt(V);
            sizeL[i] = sigma;
        }
#endif
        real max_range = max(sizeL);
        
        double scale = 200.0/max_range;
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            nslices[i] = int(sizeL[i]*scale);
        }
    }
    
    void check_table(const Array<itype> &t, bool allow_unused=false) {
        for (int j = 0; j < t.nelems; j++) {
            int v0 = t.data[j];
            if (v0 == TABLE_UNUSED) {
                if (!allow_unused) {
                    fprintf(stderr, "unfilled patch %d\n", j);
                    ASSERT2(false, "unfilled patch after dt");
                } else {
                    continue;
                }
            }
            int x = INT_TO_X(v0);
            int y = INT_TO_Y(v0);
            ASSERT2(x >= gck_xmin(wh0) && x < gck_xmax(wh0), "x out of bounds in table after dt");
            ASSERT2(y >= gck_ymin(wh0) && y < gck_ymax(wh0), "y out of bounds in table after dt");
        }

    }

    void do_after_dt() {
#if TABLE_DEBUG
        check_table(table);
#endif

        if (p->randomize_dt) {
            double T_start_randomize = wall_time();
            table_sets->compute_sets();
            for (int j = 0; j < table.nelems; j++) {
                int v0 = table.data[j];
                int xsrc = INT_TO_X(v0)-gck_xmin(wh0);
                int ysrc = INT_TO_Y(v0)-gck_ymin(wh0);
                int current = PATCH_INDEX(xsrc, ysrc);
                vector<int> *current_set = table_sets->get_set(current);
                if (current_set) {
                    int i = rand()%int(current_set->size()+1);
                    
                    if (i != 0) {
                        i--;
                        current = (*current_set)[i];
                        ASSERT2(in_bounds(current, (b0.width()-p->patch_w+1)*(b0.height()-p->patch_w+1)), "current in bounds");
                        int xcurrent = PATCH_INDEX_TO_X(current)+gck_xmin(wh0);
                        int ycurrent = PATCH_INDEX_TO_Y(current)+gck_ymin(wh0);
                        table.data[j] = XY_TO_INT(xcurrent, ycurrent);
                    }
                }
            }
            delete table_sets;
            table_sets = NULL;
            
            double T_end_randomize = wall_time();
            if (p->verbose) {
                printf("table randomize time: %f secs\n", T_end_randomize-T_start_randomize);
            }
            
            count_table("table count after randomize");
        }
    }
    
    void do_dt_raster() {
    
        /* Distance transform to fill missing values in table */
#if TABLE_VERBOSE
        printf("PatchTableFixedN distance transform\n"); fflush(stdout);
#endif

#if TABLE_OPENMP
        if (p->parallel_dt) {
            if (p->dt_threads > 0) { omp_set_num_threads(p->dt_threads); }
            else if (p->threads > 0) { omp_set_num_threads(p->threads); }
        }
        else { omp_set_num_threads(1); }
#endif
        double T_begin_dt = wall_time();
        for (int jstep = 1; jstep >= -1; jstep -= 2) {
            int jstart = 0, jend = table.nelems;
            if (jstep < 0) {
                jstart = table.nelems-1; jend = -1;
            }
            
#if TABLE_OPENMP
            #pragma omp parallel
#endif
            {
#if TABLE_VERBOSE
                printf("PatchTableFixedN distance transform allocate arrays, jstep=%d\n", jstep); fflush(stdout);
#endif
                vector<int> src_index, grid_index;
                src_index.resize(GRID_NCHANNELS);
                grid_index.resize(GRID_NCHANNELS);

#if TABLE_OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (int j0_dim0 = 0; j0_dim0 < table.sizes[0]; j0_dim0++) {
#if TABLE_VERBOSE
                    printf("PatchTableFixedN j0_dim0=%d\n", j0_dim0); fflush(stdout);
#endif
                    for (int j0 = 0; j0 < table.nelems/table.sizes[0]; j0++) {
                        int j0_full = j0_dim0 * (table.nelems/table.sizes[0]) + j0;
                        int j = jstart + j0_full*jstep;
                        
                        itype v0 = table.data[j];
                        int xbest = INT_TO_X(v0);
                        int ybest = INT_TO_Y(v0);
                        for (int i = 0; i < GRID_NCHANNELS; i++) {
                            grid_index[i] = (j/table.stride[i]) % table.sizes[i];
#if !TABLE_OPTIMIZE_DT
                            src_index[i] = grid_index[i];
#endif
                        }
                        ASSERT(table(grid_index) == v0, "expected table at calculated index to agree with original index");
                        dist_t dbest = DIST_INFINITY;
#if TABLE_OPTIMIZE_DT
                        int grid_index_int = table.index_to_int(grid_index);
#endif
                        bool dchanged = false;
                        if (v0 != TABLE_UNUSED) {
#if (!TABLE_OPTIMIZE_DT || TABLE_DT_REGULAR)
                            dbest = part->template patch_dist_to_grid<0>(wh0, xbest, ybest, grid_index);
#else
                            dbest = table_dist.data[grid_index_int];
#endif
                        }

                        if (p->dt_mode == DT_MODE_MANHATTAN || p->dt_mode == DT_MODE_SEMI_EUCLIDEAN) {
                            /* Check previous entry along each dimension in isolation (Manhattan distance transform stencil) */
                            for (int i = 0; i < GRID_NCHANNELS; i++) {
#if (!TABLE_OPTIMIZE_DT || TABLE_DT_REGULAR)
                                int orig_index = src_index[i];
                                src_index[i] -= jstep;
                                if (in_bounds(src_index[i], nslices[i])) {
                                    TABLE_TRY_DT();
                                }
                                src_index[i] = orig_index;
#else
                                int src_index_i = grid_index[i]-jstep;
                                if (in_bounds(src_index_i, nslices[i])) {
                                    int src_index_int = grid_index_int - jstep*table.stride[i];
                                    int vsrc = table.data[src_index_int];
                                    if (vsrc != TABLE_UNUSED) {
                                        int xsrc = INT_TO_X(vsrc);
                                        int ysrc = INT_TO_Y(vsrc);

                                        /* Update distance in dimension i */
                                        dist_t dcurrent = table_dist.data[src_index_int];
#if (!TABLE_OPTIMIZE_REGULAR && !TABLE_DT_REGULAR)
                                        real wh_patch = wh0(ysrc, xsrc, i);
                                        real delta_prev    = wh_patch - part->L[i]->centers[ src_index_i];
                                        dcurrent -= ((dist_t) delta_prev)*((dist_t) delta_prev);
                                        real delta_current = wh_patch - part->L[i]->centers[grid_index[i]];
                                        dcurrent += ((dist_t) delta_current)*((dist_t) delta_current);
#else
                                        dcurrent++;
#endif
                                        if (dcurrent < dbest) {
                                            dbest = dcurrent;
                                            xbest = xsrc;
                                            ybest = ysrc;
                                            dchanged = true;
                                        }
                                    }
                                }
#endif
                            }

                            if (p->dt_mode == DT_MODE_SEMI_EUCLIDEAN) {
                                /* Check previous entry along each pair of dimensions (quadratic time -- semi-Euclidean) */
                                for (int i1 = 0; i1 < GRID_NCHANNELS; i1++) {
                                    for (int i2 = 0; i2 < GRID_NCHANNELS; i2++) {
                                        int orig_index1 = src_index[i1];
                                        int orig_index2 = src_index[i2];
                                        src_index[i1] -= jstep;
                                        src_index[i2] -= jstep;
                                        if (in_bounds(src_index[i1], nslices[i1]) && in_bounds(src_index[i2], nslices[i2])) {
                                            TABLE_TRY_DT();
                                        }
                                        src_index[i1] = orig_index1;
                                        src_index[i2] = orig_index2;
                                    }
                                }
                            }
                        
                        } else if (p->dt_mode == DT_MODE_EXPTIME) {
                            /* Check previous entry along every dimension simultaneously (Euclidean distance transform stencil) */
                            for (int k = 0; k < (1<<GRID_NCHANNELS); k++) {
                                bool ok = true;
                                for (int i = 0; i < GRID_NCHANNELS; i++) {
                                    src_index[i] = grid_index[i] - (((k>>i)&1) ? jstep: 0);
                                    if (!in_bounds(src_index[i], nslices[i])) { ok = false; break; }
                                }
                                
                                if (ok) {
                                    TABLE_TRY_DT();
                                }
                            }
                        } else {
                            ASSERT2(false, "Unknown dt_mode");
                        }
#if TABLE_EXTRA_VERBOSE
                        if (dbest != DIST_INFINITY && (!in_bounds(xbest-gck_xmin(wh0), b0.width()-p->patch_w+1) || !in_bounds(ybest-gck_ymin(wh0), b0.height()-p->patch_w+1))) { fprintf(stderr, "dt source location out of bounds: %d, %d, %dx%d, patch_w=%d\n", xbest-gck_xmin(wh0), ybest-gck_ymin(wh0), b0.width()-p->patch_w+1, b0.height()-p->patch_w+1, p->patch_w); ASSERT2(false, "out of bounds"); exit(1); }
#endif

#if TABLE_OPTIMIZE_DT
                        if (dchanged) {
                            table.data[grid_index_int] = XY_TO_INT(xbest, ybest);
#if !TABLE_DT_REGULAR
                            table_dist.data[grid_index_int] = dbest;
#endif
                        }
#else
                        table(grid_index) = XY_TO_INT(xbest, ybest);
#endif
                    }
                }
            }
        }

#define PRINT_SELPATCHES() \
        { \
            vector<int> grid_index(GRID_NCHANNELS); \
            for (int i = 0; i < MIN(sel_patches.size(), 20); i++) { \
                int j = sel_patches[i]; \
                int v = table.data[j]; \
                \
                for (int k = 0; k < GRID_NCHANNELS; k++) { \
                    grid_index[k] = (j/table.stride[k]) % table.sizes[k]; \
                } \
                \
                float dist_value = part->template patch_dist_to_grid<0>(wh0, INT_TO_X(v), INT_TO_Y(v), grid_index, 0.0, true); \
                printf("patch %d: %d (%d, %d), d=%f\n", i, v, INT_TO_X(v), INT_TO_Y(v), dist_value); \
            } \
        }

#if TABLE_EXTRA_VERBOSE
        printf("PatchTableFixedN: after dt\n"); fflush(stdout);

        count_table("Count after dt");
#endif

#if (TABLE_OPTIMIZE_DT && !TABLE_DT_REGULAR)
        table_dist.resize(1);
#endif
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (raster): %f secs\n", T_end_dt-T_begin_dt);
        }
    }

#define TABLE_PROP_LABEL(i, dir) (2*(i)+1+((dir)+1)/2)

    void fill_zeros(const Array<real> &choose_wh0) {
        itype vfill = XY_TO_INT(p->patch_w-1, p->patch_w-1);
        if (allowed_patches) {
            bool found = false;
            TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER();
                vfill = XY_TO_INT(x, y);
                found = true;
                break;
            TABLE_END_FOR_ALLOWED_PATCHES();
            ASSERT2(found, "fill_zeros: found no patches that have allowed_patches set\n");
        }
        for (int i = 0; i < table.nelems; i++) {
            if (table.data[i] == TABLE_UNUSED) {
                table.data[i] = vfill;
            }
        }
    }

    void do_dt_prop() {
        double T_begin_dt = wall_time();
        
        int nmax = MAX(table.nelems - table_add_count, table_add_count);
        vector<PropElement<itype> > prop_list0(nmax);
        vector<PropElement<itype> > prop_list1(table.nelems - table_add_count);
        int current_count = 0;
        for (int j = 0; j < table.nelems; j++) {
            itype v = table.data[j];
            if (v != TABLE_UNUSED) {
                prop_list0[current_count++] = PropElement<itype>(j, v);
            }
        }
        if (current_count != table_add_count) {
            fprintf(stderr, "current_count=%d, table_add_count=%d, table.nelems=%d\n", current_count, table_add_count, table.nelems);
            ASSERT2(current_count == table_add_count, "expected current_count == table_add_count");
        }
        
        int nlabels = 1+GRID_NCHANNELS*2;
        vector<vector<PropEdge> > next_label;
        
        vector<PropEdge> next0;
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            next0.push_back(PropEdge(i, -1, nlabels));
            next0.push_back(PropEdge(i,  1, nlabels));
        }
        next_label.push_back(next0);
        
        for (int i = 0; i < GRID_NCHANNELS; i++) {
            for (int dir = -1; dir <= 1; dir += 2) {
                vector<PropEdge> next;
                next.push_back(PropEdge(i, dir, nlabels));
                for (int j = i+1; j < GRID_NCHANNELS; j++) {
                    next.push_back(PropEdge(j, -1, nlabels));
                    next.push_back(PropEdge(j,  1, nlabels));
                }
                next_label.push_back(next);
            }
        }
        ASSERT2(next_label.size() == nlabels, "expected next_label size to match nlabels");
        int prop_list_current_count = table_add_count;
        
        for (int iter = 0;; iter++) {
            if (p->dt_iters > 0 && iter > p->dt_iters) { break; }

            int even = iter%2 == 0;
            vector<PropElement<itype> > *prop_list_current = even ? &prop_list0: &prop_list1;
            vector<PropElement<itype> > *prop_list_next    = even ? &prop_list1: &prop_list0;
            if (prop_list_current_count == 0) { break; }
            int prop_list_next_count = 0;
            //prop_list_next->clear();
            //printf("iter=%d, prop_list_current size=%d, nmax=%d, table_add_count=%d, prop_list1 size=%d\n", iter, prop_list_current_count, nmax, table_add_count, table.nelems - table_add_count);
            
            for (int i = 0; i < prop_list_current_count; i++) {
                PropElement<itype> elem((*prop_list_current)[i]);
                int e_grid_index = elem.grid_index;
                ASSERT((unsigned) e_grid_index < (unsigned) table.nelems, "expected e_grid_index in bounds");
                int e_idx = elem.label;
                ASSERT((unsigned) e_idx < (unsigned) next_label.size(), "expected e_idx in bounds");
                vector<PropEdge> *e_list = &next_label[e_idx];
                
                for (int j = 0; j < (int) e_list->size(); j++) {
                    PropEdge e = (*e_list)[j];

                    int grid_index_p = e_grid_index + e.dir * table.stride[e.dim];
                    if (in_bounds(grid_index_p, table.nelems) && table.data[grid_index_p] == TABLE_UNUSED) {
                        int dim_index_p = (e_grid_index / table.stride[e.dim]) % table.sizes[e.dim];
                        dim_index_p += e.dir;
                        if (in_bounds(dim_index_p, table.sizes[e.dim])) {
                            /* Push neighboring voxel on to prop_list_next */
                            table.data[grid_index_p]       = elem.table_value; //table.data[e_grid_index];
                            (*prop_list_next)[prop_list_next_count++] = PropElement<itype>(grid_index_p, elem.table_value, e.label);
                        }
                    }
                }
            }
            
            prop_list_current_count = prop_list_next_count;
        }
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (prop): %f secs\n", T_end_dt-T_begin_dt);
        }
    }

#if TABLE_PRODUCT_QUANTIZE
    void do_dt_product_quantize() {
        double T_begin_dt = wall_time();

        if (p->product_quantize_dt_all) {
            table_add_count = (gck_xmax(wh0)-gck_xmin(wh0))*(gck_ymax(wh0)-gck_ymin(wh0));
        }
        
        Array<real> patch_arr(table_add_count, TABLE_NCHANNELS);
        vector<itype> patch_arr_index(table_add_count);
        
        int current_count = 0;
        if (p->product_quantize_dt_all) {
            for (int y = gck_ymin(wh0); y < gck_ymax(wh0); y++) {
                for (int x = gck_xmin(wh0); x < gck_xmax(wh0); x++) {
                    real *patch_row = &patch_arr(current_count, 0);
                    real *wh0_row = &wh0(y, x, 0);
                    for (int k = 0; k < TABLE_NCHANNELS; k++) {
                        patch_row[k] = wh0_row[k];
                    }
                    patch_arr_index[current_count] = XY_TO_INT(x, y);
                    
                    current_count++;
                }
            }
        } else {
            for (int j = 0; j < table.nelems; j++) {
                itype v = table.data[j];
                if (v != TABLE_UNUSED) {
                    real *patch_row = &patch_arr(current_count, 0);
                    int v_x = INT_TO_X(v);
                    int v_y = INT_TO_Y(v);
                    real *wh0_row = &wh0(v_y, v_x, 0);
                    for (int k = 0; k < TABLE_NCHANNELS; k++) {
                        patch_row[k] = wh0_row[k];
                    }
                    patch_arr_index[current_count] = v;
                    
                    current_count++;
                }
            }
        }
        if (current_count != table_add_count) {
            fprintf(stderr, "current_count=%d, table_add_count=%d, table.nelems=%d\n", current_count, table_add_count, table.nelems);
            ASSERT2(current_count == table_add_count, "expected current_count == table_add_count");
        }
        
        flann::Matrix<real> flann_patch_arr(patch_arr.data, patch_arr.height(), patch_arr.width());
//        flann::Index<flann::L2<real> > flann_patch_index(flann_patch_arr, flann::LinearIndexParams());
        flann::Index<flann::L2<real> > flann_patch_index(flann_patch_arr, flann::KDTreeSingleIndexParams());
        flann_patch_index.buildIndex();
        
        vector<real> query_vector(TABLE_NCHANNELS);
        
        vector<vector<int> > matched_indices_mat;
        vector<vector<real> > matched_dists_mat;
        flann::Matrix<real> query_matrix(&query_vector[0], 1, query_vector.size());
        matched_indices_mat.resize(1);
        matched_indices_mat[0].resize(1);
        matched_dists_mat.resize(1);
        matched_dists_mat[0].resize(1);
        flann::SearchParams search_params(p->flann_checks);
        search_params.eps = 0.0;

        double T_begin_search = wall_time();
        
        for (int j = 0; j < table.nelems; j++) {
            itype v = table.data[j];
            if (v == TABLE_UNUSED) {
                product_quantizer->get_cluster_center(j, &query_vector[0]);
                flann_patch_index.knnSearch(query_matrix, matched_indices_mat, matched_dists_mat, 1, search_params);
                int idx_kdtree = matched_indices_mat[0][0];
                ASSERT2(in_bounds(idx_kdtree, patch_arr_index.size()), "expected idx_kdtree in bounds");
                int idx = patch_arr_index[idx_kdtree];
                table.data[j] = idx;
            }
        }

        /*
        vector<PropElementProductQuantize> q(table_add_count);
        
        int current_count = 0;
#if TABLE_PRODUCT_QUANTIZE_ACTUAL_DIST
        vector<int> quantized_index(product_quantizer->product_count);
#endif
        for (int j = 0; j < table.nelems; j++) {
            itype v = table.data[j];
            if (v != TABLE_UNUSED) {
                float dist = 0;
#if TABLE_PRODUCT_QUANTIZE_ACTUAL_DIST
                int v_x = INT_TO_X(v);
                int v_y = INT_TO_Y(v);
                ASSERT2(in_bounds(v_x-gck_xmin(wh0), gck_xmax(wh0)-gck_xmin(wh0)), "expected v_x in bounds");
                ASSERT2(in_bounds(v_y-gck_ymin(wh0), gck_ymax(wh0)-gck_ymin(wh0)), "expected v_y in bounds");
                product_quantizer->quantize(&wh0(v_y, v_x, 0), &quantized_index[0]);
                for (int product_dim = 0; product_dim < product_quantizer->product_count; product_dim++) {
                    Array<real> &cluster_centers = (*product_quantizer->cluster_centers[product_dim]);
                    int dim_lo = product_quantizer->dim_lo[product_dim];
                    for (int j = 0; j < product_quantizer->dim_count[product_dim]; j++) {
                        int dim = dim_lo + j;
                        ASSERT2(in_bounds(dim, quantized_index.size()), "expected dim in bounds");
                        int icluster = quantized_index[dim];
                        ASSERT2(in_bounds(icluster, cluster_centers.height()), "expected icluster in bounds");
                        ASSERT2(in_bounds(j, cluster_centers.width()), "expected j in bounds");
                        
                        real coord2 = cluster_centers(icluster, j);
                        real coord1 = wh0(v_y, v_x, dim);
                        float delta = coord1-coord2;
                        dist += delta*delta;
                    }
                }
#endif
                q[current_count++] = PropElementProductQuantize<itype>(j, v, dist);
            }
        }
        if (current_count != table_add_count) {
            fprintf(stderr, "current_count=%d, table_add_count=%d, table.nelems=%d\n", current_count, table_add_count, table.nelems);
            ASSERT2(current_count == table_add_count, "expected current_count == table_add_count");
        }
#if TABLE_PRODUCT_QUANTIZE_ACTUAL_DIST
        make_heap(&q[0], (&q[0]) + q.size());
#endif

        while (q.size()) {
            PropElementProductQuantize e(q[0]);
            pop_heap(&q[0], &q[q.size()]);
            q.pop_back();
            ...
        }
        
        */
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (product_quantize): %f secs (%f search)\n", T_end_dt-T_begin_dt, T_end_dt-T_begin_search);
        }

    }
#endif

    void do_dt_prop1() {
        double T_begin_dt = wall_time();
        
        for (int j = 0; j < table.nelems; j++) {
            itype v = table.data[j];
            if (v != TABLE_UNUSED && !(v & TABLE_HI_MASK)) {
                for (int i = 0; i < GRID_NCHANNELS; i++) {
                    int grid_index = (j / table.stride[i]) % table.sizes[i];
                    if (grid_index > 0) {
                        int j1 = j - table.stride[i];
                        if (table.data[j1] == TABLE_UNUSED) { table.data[j1] = v | TABLE_HI_MASK; }
                    }
                    if (grid_index+1 < table.sizes[i]) {
                        int j1 = j + table.stride[i];
                        if (table.data[j1] == TABLE_UNUSED) { table.data[j1] = v | TABLE_HI_MASK; }
                    }
                }
            }
        }
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (prop1): %f secs\n", T_end_dt-T_begin_dt);
        }
    }
    
    void unmask_table() {
        for (int j = 0; j < table.nelems; j++) {
            table.data[j] = table.data[j] & TABLE_LO_MASK;
        }
    }

#if TABLE_DT_DOWNSAMPLE
    void do_dt_downsample() {
        double T_begin_dt = wall_time();
        
        vector<Array<itype> *> pyr;
        pyr.push_back(&table);
        
        vector<int> unitary(GRID_NCHANNELS, 1);
        
        while (pyr[pyr.size()-1]->sizes != unitary) {
            pyr.push_back(new Array<itype>());
            dt_downsample(*pyr[pyr.size()-2], *pyr[pyr.size()-1]);
        }
        
        for (int i = pyr.size()-1; i >= 1; i--) {
            dt_upsample(*pyr[i], *pyr[i-1]);
            delete pyr[i];
        }
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (downsample): %f secs\n", T_end_dt-T_begin_dt);
        }
    }
#endif

    void do_dt_brute() {
        double T_begin_dt = wall_time();
#if TABLE_VERBOSE
        printf("PatchTableFixedN distance transform brute\n"); fflush(stdout);
#endif

        vector<int> grid_index;
        grid_index.resize(GRID_NCHANNELS);

#if TABLE_OPENMP
        #pragma omp parallel for
#endif
        for (int j = 0; j < table.nelems; j++) {
            int xbest = gck_xmin(wh0);
            int ybest = gck_ymin(wh0);
            
            for (int i = 0; i < GRID_NCHANNELS; i++) {
                grid_index[i] = (j/table.stride[i]) % table.sizes[i];
            }
            dist_t dbest = part->template patch_dist_to_grid<0>(wh0, xbest, ybest, grid_index);

            for (int ysrc = gck_ymin(wh0); ysrc < gck_ymax(wh0); ysrc++) {
                for (int xsrc = gck_xmin(wh0); xsrc < gck_xmax(wh0); xsrc++) {
                    dist_t dcurrent = part->template patch_dist_to_grid<0>(wh0, xsrc, ysrc, grid_index);
                    if (dcurrent < dbest) {
                        dbest = dcurrent;
                        xbest = xsrc;
                        ybest = ysrc;
                    }
                }
            }

            table(grid_index) = XY_TO_INT(xbest, ybest);
        }

#if TABLE_EXTRA_VERBOSE
        printf("PatchTableFixedN: after dt brute\n"); fflush(stdout);

        count_table("Count after dt brute");
#endif

#if (TABLE_OPTIMIZE_DT && !TABLE_DT_REGULAR)
        table_dist.resize(1);
#endif
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (brute): %f secs\n", T_end_dt-T_begin_dt);
        }
    }

    void flann_knn_search(vector<real> &query_vector, vector<int> &matched_indices, vector<float> &matched_dists, int knn, flann::SearchParams &search_params) {
        flann::Matrix<real> query_matrix(&query_vector[0], 1, query_vector.size());
        static vector<vector<int> > matched_indices_mat;
        static vector<vector<float> > matched_dists_mat;
        matched_indices_mat.resize(1);
        matched_dists_mat.resize(1);
        matched_indices_mat[0].resize(knn);
        matched_dists_mat[0].resize(knn);
        ASSERT(matched_indices_mat.size() == 1, "expected 1 length for matched_indices_mat");
        ASSERT(matched_indices_mat[0].size() == knn, "expected knn length for matched_indices_mat[0]");
//        printf("knnSearch: query_matrix: %dx%d, matched_indices_mat: %dx%d, matched_dists_mat: %dx%d, knn: %d, query_vector: %d, matched_indices: %d, matched_dists: %d\n", query_matrix.rows, query_matrix.cols, matched_indices_mat.size(), matched_indices_mat[0].size(), matched_dists_mat.size(), matched_dists_mat[0].size(), knn, int(query_vector.size()), int(matched_indices.size()), int(matched_dists.size()));
        flann_index->knnSearch(query_matrix, matched_indices_mat, matched_dists_mat, knn, search_params);
        //printf("done knnSearch\n");
        for (int i = 0; i < knn; i++) {
            matched_indices[i] = matched_indices_mat[0][i];
            matched_dists[i] = matched_dists_mat[0][i];
        }
    }
    
#if TABLE_ENABLE_FLANN
    void do_dt_kdtree() {
        double T_begin_dt = wall_time();
#if TABLE_VERBOSE
        printf("PatchTableFixedN distance transform kdtree\n"); fflush(stdout);
#endif

        build_flann_index();
        
        vector<int> grid_index(TABLE_NCHANNELS);
        vector<real> query_vector(TABLE_NCHANNELS);
        vector<int> matched_indices(p->dt_knn);
        vector<float> matched_dists(p->dt_knn);
        flann::SearchParams search_params(p->flann_checks);

        int bew = gck_xmax(wh0)-gck_xmin(wh0);
        int beh = gck_ymax(wh0)-gck_ymin(wh0);

#if TABLE_OPENMP
        #pragma omp parallel for
#endif
        for (int j = 0; j < table.nelems; j++) {
            for (int i = 0; i < TABLE_NCHANNELS; i++) {
                grid_index[i] = (j/table.stride[i]) % table.sizes[i];
                query_vector[i] = part->L[i]->centers[grid_index[i]];
            }
            
            flann_knn_search(query_vector, matched_indices, matched_dists, p->dt_knn, search_params);
            int matched_index = matched_indices[rand()%p->dt_knn];
            int bx = matched_index%bew;
            int by = matched_index/bew;
            if (!in_bounds(by, beh) || !in_bounds(bx, bew)) {
                fprintf(stderr, "expected bx=%d, by=%d in range bew=%d, beh=%d", bx, by, bew, beh); ASSERT2(false, "by out of range");
            }

            table(grid_index) = XY_TO_INT(gck_xmin(wh0)+bx, gck_ymin(wh0)+by);
        }

#if TABLE_EXTRA_VERBOSE
        printf("PatchTableFixedN: after dt kdtree\n"); fflush(stdout);

        count_table("Count after dt kdtree");
#endif

#if (TABLE_OPTIMIZE_DT && !TABLE_DT_REGULAR)
        table_dist.resize(1);
#endif
        
        double T_end_dt = wall_time();
        if (p->verbose) {
            printf("table dt time (kdtree): %f secs\n", T_end_dt-T_begin_dt);
        }
    }
#endif

#define TABLE_IMG_WIDTH(wh) ((wh).width() - (p->patch_w/2) - ((p->patch_w/2)-1))
#define TABLE_IMG_HEIGHT(wh) ((wh).height() - (p->patch_w/2) - ((p->patch_w/2)-1))

    void get_padded_descriptor(const Array<in_type> &a0, Array<real> &a_image, Array<real> &a_wh, string log_message) {
        int rpad = p->patch_w/2;
        int lpad = rpad-1;
        int pad = lpad+rpad;
        int channels = 3;
        
        if (p->is_descriptor) {
            if (p->descriptor_padded) {
                a_wh.assign(a0);
            } else {
                int wh_width = a0.width()+pad+p->patch_w-1;
                int wh_height = a0.height()+pad+p->patch_w-1;
                a_wh.resize(wh_height, wh_width, a0.channels());
                a0.copy_rect(a_wh, vector<int>({0, 0, 0}), vector<int>({p->patch_w-1, p->patch_w-1, 0}), vector<int>({a0.height(), a0.width(), a0.channels()}));
            }
            a_image.resize(TABLE_IMG_HEIGHT(a_wh), TABLE_IMG_WIDTH(a_wh), channels);
        } else {
            if (p->convert_colorspace) {
                double T0 = wall_time();
                if (p->colorspace == PATCHTABLE_COLORSPACE_YUV) {
                    a0.rgb2yuv(a_image);
                } else if (p->colorspace == PATCHTABLE_COLORSPACE_LAB) {
                    a_image.assign(a0);
                    a_image.rgb2lab();
                } else {
                    fprintf(stderr, "colorspace unsupported: %d\n", p->colorspace); exit(1);
                }
                double T1 = wall_time();
                if (p->verbose) {
                    printf("%s colorspace conversion time: %f\n", log_message.c_str(), T1-T0);
                }
            } else {
                a_image.assign(a0);
            }
            const Array<in_type> &a(p->convert_colorspace ? a_image: a0);
            
            double T0_gck = wall_time();
            reduce_dim(a_image, a_wh);
            double T1_gck = wall_time();
            if (p->verbose) {
                printf("%s gck time: %f\n", log_message.c_str(), T1_gck-T0_gck);
            }
        }
#if TABLE_EXTRA_VERBOSE
        for (int y = gck_ymin(a_wh); y < gck_ymax(a_wh); y++) {
            for (int x = gck_xmin(a_wh); x < gck_xmax(a_wh); x++) {
                printf("%d, %d ", x, y);
                for (int i = 0; i < ndims; i++) {
                    printf("%.3f ", a_wh(y, x, i));
                }
                printf("\n");
            }
        }
#endif
    }
    
#if TABLE_ENABLE_ANN
    void build_ann_index() {
        if (!ann_index) {
            ann_index = build_ann_index_func<real, itype, ndims>(p, wh0, allowed_patches.get(), p->ann_bgrid, ann_index_to_pos, ann_data_points);
        }
    }
#endif

    void build_flann_index() {
#if TABLE_ENABLE_FLANN
        if (!flann_index) {
            flann_index = build_flann_index_func<real, itype, ndims>(p, wh0, allowed_patches.get(), p->flann_build_step, ann_index_to_pos);
        }
#else
        ASSERT2(false, "TABLE_ENABLE_FLANN disabled, cannot use lookup_algo LOOKUP_ALGO_KDTREE");
#endif

    }
    
    void set_nslices(const Array<real> &cluster_wh0) {
        double T0_slices = wall_time();
        set_min_max(cluster_wh0);
        if (p->nslices.size() == 0) {
            /* Determine nslices automatically */
            set_default_nslices(cluster_wh0);
        } else {
            nslices = p->nslices;
        }
        nslices = scale_slices_to_limit(nslices, p->limit);
        ASSERT2(nslices.size() == GRID_NCHANNELS, "expected nslices size to match GRID_NCHANNELS");
        if (p->verbose) {
            double T1_slices = wall_time();
            printf("table slices time: %f secs\n", T1_slices-T0_slices);
        }
    }

    PatchTableFixedN(PatchTableParams *p_, const Array<in_type> &b, Array<itype> *allowed_patches_=NULL, Array<in_type> *b0_no_copy=NULL, Array<real> *wh0_no_copy=NULL) {
#if TABLE_ENABLE_ANN
        ann_data_points = NULL;
#endif
#if TABLE_VERBOSE
        printf("PatchTableFixedN constructor (%p)\n", (void *) p_); fflush(stdout);
#endif
        if (p_->load_filename.size()) {
            load(p_, p_->load_filename);
            return;
        }
        
        bool dealloc = (b0_no_copy == NULL);
        
        double T0_construct_table = wall_time();

        table_add_count = 0;
        ASSERT2(b.dimensions() == 3, "expected 3 dimension image as input to PatchTable");
        p = p_;
        if (p->grid_ndims < 0) { p->grid_ndims = p->ndims; }
        if (p->lookup_algo == LOOKUP_ALGO_TREECANN) {
            p->limit = 1;
            p->partition_step = 10000;
            p->kcoherence = 0;
        }
        
#if TABLE_VERBOSE
        printf("copying allowed_patches\n"); fflush(stdout);
#endif
        if (allowed_patches_) {
            printf("before assign allowed_patches\n");
            allowed_patches = shared_ptr<Array<itype> >(new Array<itype>(*allowed_patches_, dealloc));
            printf("after assign allowed_patches\n");
        }
        part = NULL;
        
#if TABLE_VERBOSE
        printf("converting colorspace\n"); fflush(stdout);
#endif
        if (dealloc) {
            Array<real> b_colorspace_conv;
            Array<real> &b_image(p->is_descriptor ? b0: b_colorspace_conv);
            get_padded_descriptor(b, b_image, wh0, "table");
            if (!p->is_descriptor) {
                b0.assign(b);
            }
        } else {
            printf("before assign b0\n");
            b0.assign(*b0_no_copy, false);
            printf("after assign b0\n");
            printf("before assign wh0\n");
            wh0.assign(*wh0_no_copy, false);
            printf("after assign wh0\n");
        }
        
        int bnn_height = TABLE_IMG_HEIGHT(wh0) - p->patch_w + 1;
        int bnn_width = TABLE_IMG_WIDTH(wh0) - p->patch_w + 1;
        if (allowed_patches) {
            ASSERT2(allowed_patches->height() >= bnn_height, "expected allowed_patches height to be at least input image height-patch_w+1");
            ASSERT2(allowed_patches->width()  >= bnn_width,  "expected allowed_patches width to be at least input image width-patch_w+1");
        }

        int bew = gck_xmax(wh0)-gck_xmin(wh0);
        int beh = gck_ymax(wh0)-gck_ymin(wh0);

        /* ---------------------------------------------------------------------------------
           Calculate k-coherence
           --------------------------------------------------------------------------------- */

        if (p->kcoherence) {
#if TABLE_VERBOSE
            printf("k-coherence\n"); fflush(stdout);
#endif
            double T0_kcoherence = wall_time();

            kcoherence_set.resize(beh, bew, p->kcoherence);

            if (p->kcoherence_algo == KCOHERENCE_ALGO_FLANN || p->kcoherence_algo == KCOHERENCE_ALGO_ANN) {
#if TABLE_ENABLE_FLANN
                if (p->kcoherence_algo == KCOHERENCE_ALGO_FLANN) {
                    build_flann_index();
                }
#endif
#if TABLE_ENABLE_ANN
                if (p->kcoherence_algo == KCOHERENCE_ALGO_ANN) {
                    build_ann_index();
                }
#endif
                int knn0 = p->kcoherence+1;
                
                vector<real> query_vector(TABLE_NCHANNELS);
                vector<bool> is_valid(knn0);
                vector<int> matched_indices(knn0);
                vector<float> matched_dists(knn0);
                
#if TABLE_ENABLE_ANN
                ANNcoord ann_query[ndims];
                vector<ANNidx> ann_idx(knn0);
                vector<ANNdist> ann_dist(knn0);
#endif
#if TABLE_ENABLE_FLANN
                if (p->verbose) {
                    printf("k-coherence using flann, kcoherence_step=%d, knn=%d\n", p->kcoherence_step, p->kcoherence);
                }
                flann::SearchParams search_params(p->flann_checks);
                search_params.eps = p->flann_eps;
#endif
                
                int prev_knn = knn0;
                int min_dist2 = (p->kcoherence_min_dist*p->kcoherence_min_dist);

#if TABLE_KCOHERENCE_TRIANGLE
                PATCH_DIST_CALC_BVEC0(wh0);
                vector<pair<float, pair<int, int> > > matched_triples;
#endif
                
#if TABLE_SAVE_KCOHERENCE
                Array<float> array_kcoherence_x, array_kcoherence_y, array_kcoherence_d;
                if (p->save_kcoherence) {
                    array_kcoherence_x.resize(kcoherence_set.height(), kcoherence_set.width(), p->kcoherence);
                    array_kcoherence_y.resize(kcoherence_set.height(), kcoherence_set.width(), p->kcoherence);
                    array_kcoherence_d.resize(kcoherence_set.height(), kcoherence_set.width(), p->kcoherence);
                }
#endif
                
                TABLE_FOR_ALLOWED_PATCHES_STEP(p->kcoherence_step)
                    int y_kcoherence = (y-gck_ymin(wh0))/p->kcoherence_step;
                    int x_kcoherence = (x-gck_xmin(wh0))/p->kcoherence_step;

                    for (int i = 0; i < TABLE_NCHANNELS; i++) {
                        query_vector[i] = wh0(y, x, i);
                    }
                
                    int knn = (x/p->kcoherence_step) % p->flann_reset_step == 0 ? knn0: prev_knn;
                
                    while (1) {
                        is_valid.resize(knn);
                        matched_indices.resize(knn);
                        matched_dists.resize(knn);
                        
//                        printf("kd-tree query %d, knn=%d, kcoherence=%d at %d, %d\n", int(query_vector.size()), knn, p->kcoherence, x, y);
                        bool search_ok = false;
#if TABLE_ENABLE_FLANN
                        if (p->kcoherence_algo == KCOHERENCE_ALGO_FLANN) {
                            flann_knn_search(query_vector, matched_indices, matched_dists, knn, search_params);
                            search_ok = true;
                            for (int i = 0; i < knn; i++) {
                                if (!in_bounds(matched_indices[i], ann_index_to_pos.size())) {
                                    fprintf(stderr, "matched_indices[%d] = %d, out of bounds %d\n", i, matched_indices[i], int(ann_index_to_pos.size())); ASSERT2(false, "out of bounds");
                                }
                                matched_indices[i] = ann_index_to_pos[matched_indices[i]];
                            }
                        }
#endif
#if TABLE_ENABLE_ANN
                        if (p->kcoherence_algo == KCOHERENCE_ALGO_ANN) {
                            ann_idx.resize(knn);
                            ann_dist.resize(knn);
                            for (int i = 0; i < TABLE_NCHANNELS; i++) {
                                ann_query[i] = query_vector[i];
                            }
                            ann_index->annkSearch(ann_query, knn, &ann_idx[0], &ann_dist[0], p->ann_eps);
                            for (int i = 0; i < knn; i++) {
                                if (!in_bounds(ann_idx[i], ann_index_to_pos.size())) {
                                    fprintf(stderr, "ann_idx[%d] = %d, out of bounds %d\n", i, ann_idx[i], ann_index_to_pos.size()); ASSERT2(false, "out of bounds");
                                }
                                matched_indices[i] = ann_index_to_pos[ann_idx[i]];
                            }
                            search_ok = true;
                        }
#endif
                        ASSERT2(search_ok, "ann/flann search failed. Enable TABLE_ENABLE_FLANN or TABLE_ENABLE_ANN");
                        
                        int nvalid = 0;
                        for (int i = 0; i < knn; i++) {
                            int matched_index = matched_indices[i];
                            int bx = INT_TO_X(matched_index), by = INT_TO_Y(matched_index);
                            if (!in_bounds(by, beh) || !in_bounds(bx, bew)) {
                                fprintf(stderr, "expected bx=%d, by=%d in range bew=%d, beh=%d, bgrid=%d, i=%d/%d, matched_indices[i]=%d\n", bx, by, bew, beh, p->ann_bgrid, i, knn, matched_indices[i]); ASSERT2(false, "by out of range");
                            }
                            int dx = bx+gck_xmin(wh0) - x;
                            int dy = by+gck_ymin(wh0) - y;
                            int d = dx*dx+dy*dy;
                            is_valid[i] = d >= min_dist2;
                            
                            if (p->incoherence) {
                                int dx_L[] = { -1,  0,  0 };
                                int dy_L[] = {  0, -1,  0 };
                                int j_delta_max = 2;
                                if (p->incoherence == 2) { j_delta_max = 3; }

                                for (int j_delta = 0; j_delta < j_delta_max; j_delta++) {
                                    int xp = x_kcoherence + dx_L[j_delta];
                                    int yp = y_kcoherence + dy_L[j_delta];
                                    if (xp >= 0 && yp >= 0) {
                                        for (int k = 0; k < (j_delta == 2 ? i: knn); k++) {
                                            int bx_p = 0, by_p = 0;
                                            if (j_delta < 2) {
                                                int v_p = TABLE_KCOHERENCE_POS(kcoherence_set(yp, xp, k));
                                                bx_p = INT_TO_X(v_p);
                                                by_p = INT_TO_Y(v_p);
                                            } else {
                                                if (!is_valid[k]) { continue; }
                                                bx_p = INT_TO_X(matched_indices[k]);
                                                by_p = INT_TO_Y(matched_indices[k]);
                                            }
                                            bx_p -= dx_L[j_delta] * p->kcoherence_step;
                                            by_p -= dy_L[j_delta] * p->kcoherence_step;
                                            
                                            int dx_p = bx - bx_p;
                                            int dy_p = by - by_p;
                                            int d_p = dx_p*dx_p + dy_p*dy_p;
                                            if (d_p < p->incoherence_min_dist*p->incoherence_min_dist) {
                                                is_valid[i] = false;
                                                break;
                                            }
                                        }
                                    }
                                    if (!is_valid[i]) { break; }
                                }
                            }
                                                                
                            if (is_valid[i]) {
                                nvalid++;
                                if (nvalid >= p->kcoherence) { break; }
                            }
                        }
                    
                        if (nvalid < p->kcoherence) {
                            knn *= 2;
                        } else {
#if TABLE_KCOHERENCE_TRIANGLE
                            /* Sort matched patches by exact distance */
                            TABLE_PREPARE_B_DIST();
                            for (int i = 0; i < knn; i++) {
                                int v = matched_indices[i];
                                int xsrc = INT_TO_X(v);
                                int ysrc = INT_TO_Y(v);
                                matched_dists[i] = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                            }
                            
                            matched_triples.resize(matched_indices.size());
                            for (int i = 0; i < knn; i++) {
                                matched_triples[i] = std::make_pair(matched_dists[i], std::make_pair(matched_indices[i], is_valid[i]));
                            }
                            sort(matched_triples.begin(), matched_triples.end());
                            for (int i = 0; i < knn; i++) {
                                matched_dists[i] = matched_triples[i].first;
                                matched_indices[i] = matched_triples[i].second.first;
                                is_valid[i] = matched_triples[i].second.second;
                            }
#endif
                            
                            int current_valid = 0;
                            for (int i = 0; i < knn; i++) {
                                int matched_index = matched_indices[i];
                                int bx = INT_TO_X(matched_index);
                                int by = INT_TO_Y(matched_index);
                                if (is_valid[i]) {
                                    ASSERT(bx != x-gck_xmin(wh0) || by != y-gck_ymin(wh0), "expected k-coherence match differing from source coord");
                                    kcoherence_type *k_p = &kcoherence_set(y_kcoherence, x_kcoherence, current_valid);
                                    TABLE_KCOHERENCE_POS(*k_p) = XY_TO_INT(bx, by);
#if TABLE_KCOHERENCE_TRIANGLE
                                    TABLE_KCOHERENCE_DIST(*k_p) = matched_dists[i];
#endif
#if TABLE_SAVE_KCOHERENCE
                                    if (p->save_kcoherence) {
                                        array_kcoherence_x(y_kcoherence, x_kcoherence, current_valid) = bx;
                                        array_kcoherence_y(y_kcoherence, x_kcoherence, current_valid) = by;
                                        array_kcoherence_d(y_kcoherence, x_kcoherence, current_valid) = matched_dists[i];
                                    }
#endif
#if TABLE_VERBOSE
                                    printf("k-coherence set %d %d %d: %d %d\n", y_kcoherence, x_kcoherence, current_valid, bx, by);
#endif
                                    current_valid++;
                                    if (current_valid >= p->kcoherence) { break; }
                                }
                            }
                            ASSERT2(current_valid == p->kcoherence, "expected current_valid to equal kcoherence");
                            break;
                        }
                    }
                    prev_knn = knn;
                
                TABLE_END_FOR_ALLOWED_PATCHES()

#if TABLE_SAVE_KCOHERENCE
                if (p->save_kcoherence) {
                    save_color_image<float>(array_kcoherence_x, "kcoherence_x.pfm");
                    save_color_image<float>(array_kcoherence_y, "kcoherence_y.pfm");
                    save_color_image<float>(array_kcoherence_d, "kcoherence_d.pfm");
                }
#endif

#if TABLE_DEBUG
                check_kcoherence();
#endif
            } else if (p->kcoherence_algo == KCOHERENCE_ALGO_TREECANN) {
                Array<double> kcoherence_nnf;
                p->pm_knn = p->kcoherence;
                p->pm_min_dist = p->kcoherence_min_dist;
                
#if TABLE_VERBOSE
                printf("calling treecann\n"); fflush(stdout);
                printf("allowed_patches=%p\n", allowed_patches.get()); fflush(stdout);
#endif
                if (p->is_descriptor) {
                    TreeCANN<real, 1, ndims, 1, itype>(p, b0, wh0, allowed_patches.get(), p->kcoherence).lookup(b0, wh0, kcoherence_nnf, allowed_patches.get());
                } else {
                    TreeCANN<real, 0, ndims, 1, itype>(p, b0, wh0, allowed_patches.get(), p->kcoherence).lookup(b0, wh0, kcoherence_nnf, allowed_patches.get());
                }
                
                if (kcoherence_nnf.width() != kcoherence_set.width() || kcoherence_nnf.height() != kcoherence_set.height() || kcoherence_nnf.sizes[2] != p->kcoherence || kcoherence_set.sizes[2] != p->kcoherence) {
                    fprintf(stderr, "kcoherence_nnf sizes %s do not match kcoherence_set sizes %s\n", vector_to_str_int(kcoherence_nnf.sizes).c_str(),
                            vector_to_str_int(kcoherence_set.sizes).c_str());
                    ASSERT2(false, "size mismatch");
                }
                for (int y = 0; y < beh; y++) {
                    for (int x = 0; x < bew; x++) {
                        if (allowed_patches && ((*allowed_patches)(y, x)) != p->allowed_index) { continue; }
                        for (int k = 0; k < p->kcoherence; k++) {
                            int k_bx = kcoherence_nnf(y, x, k, NNF_X);
                            int k_by = kcoherence_nnf(y, x, k, NNF_Y);
#if TABLE_DEBUG
                            if (!in_bounds(k_bx, bew) || !in_bounds(k_by, beh)) {
                                fprintf(stderr, "kcoherence NN out of bounds: %d, %d, %d => %dx%d\n", x, y, k, k_bx, k_by);
                                ASSERT2(false, "out of bounds");
                            }
                            if (allowed_patches && ((*allowed_patches)(k_by, k_bx)) != p->allowed_index) {
                                fprintf(stderr, "kcoherence NN points to a disallowed patch: %d, %d, %d => %d, %d\n", x, y, k, k_bx, k_by);
                            }
#endif
                            TABLE_KCOHERENCE_POS(kcoherence_set(y, x, k)) = XY_TO_INT(k_bx, k_by);
#if TABLE_KCOHERENCE_TRIANGLE
                            TABLE_KCOHERENCE_DIST(kcoherence_set(y, x, k)) = kcoherence_nnf(y, x, k, NNF_DIST);
#endif
                        }
                    }
                }
                
            } else if (p->kcoherence_algo == KCOHERENCE_ALGO_PM) {
                Array<double> kcoherence_nnf;
                p->pm_knn = p->kcoherence;
                p->pm_min_dist = p->kcoherence_min_dist;
                
#if TABLE_VERBOSE
                printf("calling patchmatch\n"); fflush(stdout);
                printf("allowed_patches=%p\n", allowed_patches.get()); fflush(stdout);
#endif
                if (p->is_descriptor) {
                    patchmatch<real, 1, TABLE_NCHANNELS, 1>(p, b0, b0, wh0, wh0, kcoherence_nnf, allowed_patches.get(), allowed_patches.get(), p->kcoherence_step);
                } else {
                    patchmatch<real, 0, TABLE_NCHANNELS, 1>(p, b0, b0, wh0, wh0, kcoherence_nnf, allowed_patches.get(), allowed_patches.get(), p->kcoherence_step);
                }

                if (kcoherence_nnf.width() != kcoherence_set.width() || kcoherence_nnf.height() != kcoherence_set.height() || kcoherence_nnf.sizes[2] != p->kcoherence || kcoherence_set.sizes[2] != p->kcoherence) {
                    fprintf(stderr, "kcoherence_nnf sizes %s do not match kcoherence_set sizes %s\n", vector_to_str_int(kcoherence_nnf.sizes).c_str(),
                                                                                                      vector_to_str_int(kcoherence_set.sizes).c_str());
                    ASSERT2(false, "size mismatch");
                }
                for (int y = 0; y < beh; y += p->kcoherence_step) {
                    for (int x = 0; x < bew; x += p->kcoherence_step) {
                        int y_lo = y/p->kcoherence_step;
                        int x_lo = x/p->kcoherence_step;
                        if (allowed_patches && ((*allowed_patches)(y, x)) != p->allowed_index) { continue; }
                        for (int k = 0; k < p->kcoherence; k++) {
                            int k_bx = kcoherence_nnf(y_lo, x_lo, k, NNF_X);
                            int k_by = kcoherence_nnf(y_lo, x_lo, k, NNF_Y);
#if TABLE_DEBUG
                            if (!in_bounds(k_bx, bew) || !in_bounds(k_by, beh)) {
                                fprintf(stderr, "kcoherence NN out of bounds: %d, %d, %d => %dx%d\n", x, y, k, k_bx, k_by);
                                ASSERT2(false, "out of bounds");
                            }
                            if (allowed_patches && ((*allowed_patches)(k_by, k_bx)) != p->allowed_index) {
                                fprintf(stderr, "kcoherence NN points to a disallowed patch: %d, %d, %d => %d, %d\n", x, y, k, k_bx, k_by);
                            }
#endif
                            TABLE_KCOHERENCE_POS(kcoherence_set(y_lo, x_lo, k)) = XY_TO_INT(k_bx, k_by);
                            ASSERT2(!TABLE_KCOHERENCE_TRIANGLE, "TABLE_KCOHERENCE_TRIANGLE not supported for kcoherence_algo pm");
                        }
                    }
                }
            } else {
                ASSERT2(false, "unknown kcoherence_algo");
            }
            /*
            for (int y = 0; y < kcoherence_set.height(); y++) {
                for (int x = 0; x < kcoherence_set.width(); x++) {
                    for (int k = 0; k < p->kcoherence; k++) {
                        int v = kcoherence_set(y, x, k);
                        int bx = INT_TO_X(v);
                        int by = INT_TO_Y(v);
                        if (!in_bounds(bx, bew) || !in_bounds(by, beh)) {
                            fprintf(stderr, "kcoherence out of bounds: %d, %d, b: %d, %d, bounds %dx%d\n", x, y, bx, by, bew, beh); ASSERT2(false, "out of bounds");
                        }
                    }
                }
            }
            */
            if (p->verbose) {
                printf("table kcoherence time: %f\n", wall_time()-T0_kcoherence);
            }
        }
        
        /* ---------------------------------------------------------------------------------
           Compute k-means clusters
           --------------------------------------------------------------------------------- */

        Array<real> cluster_wh0;
#if TABLE_CLUSTER_KMEANS
        const Array<real> &choose_wh0(p->cluster_kmeans ? cluster_wh0: wh0);
        
        cv::Mat cluster_labels, cluster_centers;
        
        if (p->cluster_kmeans) {
            double T0_cluster_kmeans = wall_time();
            
            int point_count = 0;
            TABLE_FOR_ALLOWED_PATCHES() {
                point_count++;
            } TABLE_END_FOR_ALLOWED_PATCHES();

            cv::Mat cluster_points(point_count, GRID_NCHANNELS, CV_32FC1);
            
            int current_point = 0;
            TABLE_FOR_ALLOWED_PATCHES() {
                float *row = cluster_points.ptr<float>(current_point);
                real *wh0_p = &wh0.get_nearest(y, x, 0);
                for (int j = 0; j < GRID_NCHANNELS; j++) {
                    row[j] = wh0_p[j];
                }
                current_point++;
            } TABLE_END_FOR_ALLOWED_PATCHES();
            ASSERT2(current_point == point_count, "expected current_point == point_count");
            
            cv::TermCriteria criteria(cv::TermCriteria::EPS+cv::TermCriteria::COUNT, p->kmeans_max_iters, p->kmeans_eps);
            double T0_kmeans = wall_time();
            kmeans(cluster_points, p->cluster_count, cluster_labels, criteria, p->kmeans_attempts, cv::KMEANS_PP_CENTERS, cluster_centers);
            double T_kmeans = wall_time() - T0_kmeans;
            ASSERT2(cluster_labels.rows == point_count, "expected labels.rows == point_count");
            ASSERT2(cluster_labels.cols == 1, "expected labels.cols == 1");
            ASSERT2(cluster_centers.rows == p->cluster_count, "expected clusters.rows == cluster_count");
            ASSERT2(cluster_centers.cols == GRID_NCHANNELS, "expected centers.cols == GRID_NCHANNELS");
            
            // Want gck_ymax() = gck_ymin() + 1
            //      (h - pw+1) = (pw-1) + 1
            //      h = 2*pw-1
            cluster_wh0.resize(2*p->patch_w-1, 2*p->patch_w-2 + p->cluster_count, wh0.channels());
            sub_upper_left = XY_TO_INT(gck_xmin(cluster_wh0), gck_ymin(cluster_wh0));

            ASSERT2(gck_ymax(cluster_wh0) == gck_ymin(cluster_wh0)+1, "expected gck_ymax() == gck_ymin()+1");
            for (int y = 0; y < 1; y++) {
                for (int x = 0; x < p->cluster_count; x++) {
                    float *centers_row = cluster_centers.ptr<float>(x);
                    for (int i = 0; i < wh0.channels(); i++) {
                        cluster_wh0(gck_ymin(cluster_wh0), gck_xmin(cluster_wh0)+x, i) = centers_row[i];
                    }
                }
            }
            printf("table cluster kmeans time: %f\n", wall_time()-T0_cluster_kmeans);
        }
#else
        const Array<real> &choose_wh0(wh0);
#endif

        /* ---------------------------------------------------------------------------------
           Build table 
           --------------------------------------------------------------------------------- */
        
        bool table_built = false;

#if TABLE_PRODUCT_QUANTIZE
        if (p->product_quantize) {
            /* Build product quantized table */
            table_built = true;
            set_nslices(choose_wh0);
            
            product_quantizer = make_shared<ProductQuantizer<real, itype> >(p, wh0, nslices, allowed_patches.get());
            
            /* Populate table with patches */
            table.resize(nslices);
            table.clear(TABLE_UNUSED);
            if (!p->product_quantize_dt_all) {
                TABLE_FOR_ALLOWED_PATCHES() {
                    int patch_index = product_quantizer->quantize(&wh0(y, x, 0));   /* TODO: Optimize this method */
                    if (!in_bounds(patch_index, table.nelems)) {
                        fprintf(stderr, "patch_index (%d) not in table bounds %d\n", patch_index, table.nelems); ASSERT2(false, "out of bounds");
                    }
                    
                    int v = table.data[patch_index];
                    if (v == TABLE_UNUSED) {
                        table.data[patch_index] = XY_TO_INT(x, y);
                        if (TABLE_DT_REGULAR) {
                            table_add_count++;
                        }
                    }
                } TABLE_END_FOR_ALLOWED_PATCHES();
            }
            
            do_dt_product_quantize();
            check_table(table);
        }
#endif
        if (!table_built && (p->lookup_algo != LOOKUP_ALGO_KDTREE || p->kdtree_add_unique)) {
            /* Build ordinary table */
            bool unique = p->populate_nearest && p->randomize_dt;
            if (p->randomize_dt) {
                table_sets = new AdjacencySet((bnn_width)*(bnn_height), unique);
            }
#if TABLE_VERBOSE
            printf("PatchTableFixedN reduce_dim\n"); fflush(stdout);
#endif

            set_nslices(choose_wh0);
            
#if TABLE_VERBOSE
            printf("PatchTableFixedN constructing partitions\n"); fflush(stdout);
#endif

            double T_begin_part = wall_time();
            part = new PatchPartition<real, itype>(p, choose_wh0, nslices, false, allowed_patches.get(), min_val, max_val);
            double T_end_part = wall_time();
            if (p->verbose) {
                printf("table part time: %f secs\n", T_end_part-T_begin_part); fflush(stdout);
            }
            
            vector<int> orig_patch_index;
            orig_patch_index.resize(nslices.size());
            
            if (p->verbose) {
                for (int i = 0; i < (int) GRID_NCHANNELS; i++) {
                    printf("%d: [%f, %f], range=%f, nslices=%d\n", i, double(part->L[i]->min_val), double(part->L[i]->max_val), double(part->L[i]->max_val-part->L[i]->min_val), nslices[i]);
                }
                fflush(stdout);
            }

            /* Populate table with patches */
            if (p->lookup_algo != LOOKUP_ALGO_TREECANN) {
                double T_begin_table = wall_time();
                table.resize(nslices);
                table.clear(TABLE_UNUSED);
#if TABLE_POPULATE_FAST
#define TABLE_POPULATE_INNER_BLOCK() \
                    int patch_index = 0; \
                    for (int i = 0; i < GRID_NCHANNELS; i++) { \
                        patch_index += part->L[i]->get_bin(choose_wh0(y, x, i)) * table.stride[i]; \
                    } \
                    \
                    ASSERT(in_bounds(patch_index, table.nelems), "patch_index not in bounds"); \
                    int v = table.data[patch_index]; \
                    if (v == TABLE_UNUSED) { \
                        table.data[patch_index] = XY_TO_INT(x, y); \
                        if (TABLE_DT_REGULAR) { \
                            table_add_count++; \
                        } \
                    }

                if (p->populate_random) {
                    vector<int> patch_order(gck_ew(wh0)*gck_eh(wh0));
                    int patch_order_count = 0;
                    TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER();
                        patch_order[patch_order_count++] = XY_TO_INT(x, y);
                    TABLE_END_FOR_ALLOWED_PATCHES();
                    std::random_shuffle(&patch_order[0], (&patch_order[0])+patch_order_count);
                    
                    for (int i = 0; i < patch_order_count; i++) {
                        int pos = patch_order[i];
                        int x = INT_TO_X(pos), y = INT_TO_Y(pos);
                        TABLE_POPULATE_INNER_BLOCK();
                    }
                } else {
                    TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER();
                        TABLE_POPULATE_INNER_BLOCK();
                    TABLE_END_FOR_ALLOWED_PATCHES();
                }
#else       /* Not TABLE_POPULATE_FAST */
                vector<int> lo, hi, grid_index;
                lo.resize(GRID_NCHANNELS);
                hi.resize(GRID_NCHANNELS);
                grid_index.resize(GRID_NCHANNELS);
                
#if (TABLE_OPTIMIZE_DT && !TABLE_DT_REGULAR)
                table_dist.resize(nslices);
#endif
                TABLE_FOR_ALLOWED_PATCHES_OPTIONAL_CLUSTER()
                
                    part->get_patch_index(choose_wh0, x, y, orig_patch_index);
#if TABLE_VERBOSE
                    if (x % 10 == 0 && y % 10 == 0) {
                        printf("%d, %d => ", x, y);
                        for (int i = 0; i < (int) wh0.channels(); i++) {
                            printf("%f, ", double(wh0(y, x, i)));
                        }
                        printf(" => ");
                        for (int i = 0; i < (int) orig_patch_index.size(); i++) {
                            printf("%d, ", orig_patch_index[i]);
                        }
                        printf("\n");
                        fflush(stdout);
                    }
#endif
                    if (!p->populate_nearest) {
                        for (int i = 0; i < (int) orig_patch_index.size(); i++) {
                            lo[i] = orig_patch_index[i];
                            hi[i] = lo[i]+1;
                            if (hi[i] >= nslices[i]) { hi[i] = nslices[i]-1; lo[i] = hi[i]-1; }
                            if (lo[i] < 0) { lo[i] = 0; }
                            grid_index[i] = orig_patch_index[i];
                        }
                    }
                    
                    /* Add patch to adjacent voxels */
                    if (p->populate_nearest) {
                        add_patch(x, y, orig_patch_index, p, unique);
                    } else if (!p->populate_exptime) {                               /* Add +- neighbor along each dimension in isolation (linear time) */
                            for (int i = 0; i < GRID_NCHANNELS; i++) {
                                int orig_index = grid_index[i];
                                grid_index[i] = lo[i];
                                add_patch(x, y, grid_index, p, unique);
                                grid_index[i] = hi[i];
                                add_patch(x, y, grid_index, p, unique);
                                grid_index[i] = orig_index;
                            }
                    } else {                                                         /* Add +- neighbors along all dimensions jointly (exp time) */
                        for (int j = 0; j < (1<<int(GRID_NCHANNELS)); j++) {
                            for (int i = 0; i < GRID_NCHANNELS; i++) {
                                grid_index[i] = ((j>>i)&1) ? hi[i]: lo[i];
                            }
                            add_patch(x, y, grid_index, p, unique);
                        }
                    }
                TABLE_END_FOR_ALLOWED_PATCHES()
#endif          /* not TABLE_POPULATE_FAST */

                double T_end_table = wall_time();
                if (p->verbose) {
#if TABLE_DT_REGULAR
                    printf("table populate time: %f secs (%d added, %dx%d is effective size)\n", T_end_table-T_begin_table, table_add_count, gck_xmax(wh0)-gck_xmin(wh0), gck_ymax(wh0)-gck_ymin(wh0)); fflush(stdout);
#else
                    printf("table populate time: %f secs\n", T_end_table-T_begin_table); fflush(stdout);
#endif
                }
            }

#if TABLE_EXTRA_VERBOSE
            printf("table before allowed_patches\n"); fflush(stdout);
#endif

#if TABLE_EXTRA_VERBOSE
            printf("table before count\n"); fflush(stdout);
#endif

            count_table();
#if TABLE_EXTRA_VERBOSE
            printf("table before sel_patches\n"); fflush(stdout);
            vector<int> sel_patches;
            for (int j = 0; j < table.nelems; j++) {
                if (table.data[j] != TABLE_UNUSED) {
                    sel_patches.push_back(j);
                }
            }

            printf("table before print_selpatches\n"); fflush(stdout);

            PRINT_SELPATCHES();
            printf("table before second count\n"); fflush(stdout);
            count_table("Second count");
#endif

            if (p->lookup_algo != LOOKUP_ALGO_TREECANN) {
                bool zeros_filled = false;
                if (p->run_dt) {
                    if (p->dt_algo == DT_ALGO_RASTER) {
                        do_dt_raster();
                        zeros_filled = true;
                    } else if (p->dt_algo == DT_ALGO_PROP) {
                        do_dt_prop();
                        zeros_filled = (p->dt_iters <= 0);
                    } else if (p->dt_algo == DT_ALGO_BRUTE) {
                        do_dt_brute();
                        zeros_filled = true;
                    } else if (p->dt_algo == DT_ALGO_KDTREE) {
#if TABLE_ENABLE_FLANN
                        do_dt_kdtree();
#else
                        fprintf(stderr, "TABLE_ENABLE_FLANN was set to 0 when compiling, cannot use kdtree");
#endif
                        zeros_filled = true;
                    } else if (p->dt_algo == DT_ALGO_DOWNSAMPLE) {
                        do_dt_downsample();
                        zeros_filled = true;
                    } else if (p->dt_algo == DT_ALGO_HYBRID) {
                        do_dt_prop1();
                        do_dt_downsample();
                        unmask_table();
                        zeros_filled = true;
					}
					else if (p->dt_algo == DT_ALGO_EUCLIDEAN){
						do_dt_euclidean();
						zeros_filled = true;
					}else {
                        fprintf(stderr, "invalid dt_algo\n"); exit(1);
                    }
                }
                if (!zeros_filled) {
                    double T0_fill_zeros = wall_time();
                    fill_zeros(choose_wh0);
                    printf("table fill_zeros: %f secs\n", wall_time()-T0_fill_zeros);
                }
                do_after_dt();
            }
            
    /*        if (TABLE_RANDOMIZE_LOOKUP) {
                delete part;
                part = new PatchPartition<real, itype>(p, wh0, nslices, true);
            }*/
        }

#if TABLE_OPENMP
        if (p->threads > 0) { omp_set_num_threads(p->threads); }
#endif

#if TABLE_ENABLE_ANN
        if (p->ann) {
            build_ann_index();
        }
#endif
        
        if (p->lookup_algo == LOOKUP_ALGO_KDTREE) {
            build_flann_index();
        }
        
        if (p->lookup_algo == LOOKUP_ALGO_TREECANN) {
            treecann = make_shared<TreeCANN<real, 1, ndims, 0, itype> >(p, b0, wh0, allowed_patches.get());
        }
        
        if (p->verbose) {
            printf("table total precomputation time: %f secs\n", wall_time()-T0_construct_table);
            printf("\n");
        }

        if (p->save_filename.size()) {
            save(p->save_filename);
        }

#if TABLE_CLUSTER_KMEANS
        if (p->cluster_kmeans) {
            sub_allowed_patches.resize(bnn_height, bnn_width);
            sub_allowed_patches.clear(-1);
            
            int label_count = 0;
            TABLE_FOR_ALLOWED_PATCHES() {
                int lbl = cluster_labels.at<int>(label_count);
                sub_allowed_patches(y-gck_ymin(wh0), x-gck_xmin(wh0)) = lbl;
                label_count++;
            } TABLE_END_FOR_ALLOWED_PATCHES();

            for (int i = 0; i < p->cluster_count; i++) {
                if (p->verbose) {
                    printf("====================================================================\n");
                    printf("Building sub-table %d/%d\n", i, p->cluster_count);
                    printf("====================================================================\n");
                    printf("\n");
                }
                shared_ptr<PatchTableParams> pcopy = make_shared<PatchTableParams>(*p);
                pcopy->cluster_kmeans = false;
                pcopy->kcoherence = 0;
                pcopy->allowed_index = i;
                pcopy->limit = p->cluster_limit;
                sub_tables.push_back(make_shared<PatchTableFixedN<real, in_type, ndims, itype> >(pcopy.get(), b, &sub_allowed_patches, &b0, &wh0));
            }
        }
#endif

#if TABLE_MULTI_TABLES
        if (p->ntables > 1) {
            for (int i = 1; i < p->ntables; i++) {
                tables_params.push_back(make_shared<PatchTableParams>(*p));
                PatchTableParams *pcopy = tables_params[tables_params.size()-1].get();
                double mul_factor = 1.0 / pow(p->table_ratio, i);
                pcopy->limit = int(pcopy->limit*mul_factor);
                pcopy->kcoherence = 0;
                pcopy->ntables = 1;
                
                tables.push_back(make_shared<PatchTableFixedN<real, in_type, ndims, itype> >(pcopy, b, allowed_patches_));
            }
        }
#endif
    }
    
    ~PatchTableFixedN() {
        if (p->verbose) {
            printf("~PatchTableFixedN\n");
        }
        delete part;
        if (p->verbose) {
            printf("Finished ~PatchTableFixedN\n");
        }
#if TABLE_ENABLE_ANN
        if (ann_data_points) {
            annDeallocPts(ann_data_points);
            ann_data_points = NULL;
        }
#endif
    }

#define CHECK_NNF_BOUNDS(x, y, xsrc, ysrc) \
    if (TABLE_DEBUG && !in_bounds(ysrc, b0.height()-p->patch_w+1)) { fprintf(stderr, "x, y=%d, %d, ysrc %d out of bounds, xsrc=%d, h=%d, w=%d, patch_w=%d, eh=%d, ew=%d, line %d\n", x, y, ysrc, xsrc, b0.height(), b0.width(), p->patch_w, b0.height()-p->patch_w+1, b0.width()-p->patch_w+1, __LINE__); exit(1); } \
    if (TABLE_DEBUG && !in_bounds(xsrc, b0.width()-p->patch_w+1)) { fprintf(stderr, "x, y=%d, %d, xsrc %d out of bounds, ysrc=%d, h=%d, w=%d, patch_w=%d, eh=%d, ew=%d, line %d\n", x, y, xsrc, ysrc, b0.height(), b0.width(), p->patch_w, b0.height()-p->patch_w+1, b0.width()-p->patch_w+1, __LINE__); exit(1); } \
    if (TABLE_DEBUG && allowed_patches && (*allowed_patches)(ysrc, xsrc) != p->allowed_index) { fprintf(stderr, "x, y=%d, %d, xsrc=%d, %d, not allowed_patch, line %d\n", x, y, xsrc, ysrc, __LINE__); exit(1); }

#if !TABLE_PROP_EXACT
    void lookup_dist(const Array<real> &a, Array<real> &wh, Array<double> &ann, bool calc_exact_dist=true) {
        int ann_h = TABLE_IMG_HEIGHT(wh)-p->patch_w+1;
        int ann_w = TABLE_IMG_WIDTH(wh)-p->patch_w+1;

        if (p->is_descriptor) { calc_exact_dist = false; }
        
        double T4 = wall_time();
        PATCH_DIST_CALC_BVEC0(wh0);
#if TABLE_OPENMP
#pragma omp parallel for
#endif
        for (int y = 0; y < ann_h; y++) {
            for (int x = 0; x < ann_w; x++) {
                int xsrc = ann(y, x, NNF_X);
                int ysrc = ann(y, x, NNF_Y);
                CHECK_NNF_BOUNDS(x, y, xsrc, ysrc);
                
                PATCH_DIST_CALC_AVEC(wh);

                real dsrc = 0;
                if (calc_exact_dist) {
                    for (int dy = 0; dy < p->patch_w; dy++) {
                        for (int dx = 0; dx < p->patch_w; dx++) {
                            for (int c = 0; c < 3; c++) {
                                real delta = b0(ysrc+dy, xsrc+dx, c) - a(y+dy, x+dx, c);
                                dsrc += delta*delta;
                            }
                        }
                    }
                } else {
                    dsrc = patch_dist_approx<real, TABLE_NCHANNELS>(wh0, avec, bvec0, xsrc, ysrc);
                }
                ann(y, x, NNF_DIST) = dsrc;
            }
        }
        double T5 = wall_time();
        if (p->verbose) {
            printf("    lookup dist time: %f (calc_exact_dist=%d)\n", T5-T4, calc_exact_dist);
        }
    }
#endif

    void check_nnf(const Array<double> &ann, string s="checking bounds") {
        printf("check_nnf: %s\n", s.c_str());
        int ann_h = ann.height();
        int ann_w = ann.width();
        for (int y = 0; y < ann_h; y++) {
            for (int x = 0; x < ann_w; x++) {
                if (!in_bounds(ann(y, x, NNF_X), b0.width()-p->patch_w+1)) { fprintf(stderr, "p=%p, ann xsrc out of bounds at %d, %d => %d, %d (%d x %d, reduced %d x %d, patch_w=%d, p=%p)\n", (void *) p, x, y, int(ann(y, x, NNF_X)), int(ann(y, x, NNF_Y)), b0.width(), b0.height(), b0.width()-p->patch_w+1, b0.height()-p->patch_w+1, p->patch_w, (void *) p); exit(1); }
                if (!in_bounds(ann(y, x, NNF_Y), b0.height()-p->patch_w+1)) { fprintf(stderr, "ann ysrc out of bounds at %d, %d => %d, %d (%d x %d, reduced %d x %d, patch_w=%d, p=%p)\n", x, y, int(ann(y, x, NNF_X)), int(ann(y, x, NNF_Y)), b0.width(), b0.height(), b0.width()-p->patch_w+1, b0.height()-p->patch_w+1, p->patch_w, (void *) p); exit(1); }
            }
        }
    }

#define TABLE_COPY_PREV_NNF() \
    double *ann_p = &ann_row[x*3]; \
    double *ann_prev_p = &ann_prev_row[x*3]; \
    \
    int xsrc = ann_prev_p[NNF_X]; \
    int ysrc = ann_prev_p[NNF_Y]; \
    if (p->sanitize_input) { \
        if (xsrc < 0) { xsrc = 0; } \
        else if (xsrc >= b0_ew) { xsrc = b0_ew-1; } \
        if (ysrc < 0) { ysrc = 0; } \
        else if (ysrc >= b0_eh) { ysrc = b0_eh-1; } \
    } \
    ASSERT(in_bounds(xsrc, b0_ew) && in_bounds(ysrc, b0_eh), "expected xsrc, ysrc in bounds in prev_nnf mode"); \
    ann_p[NNF_X] = xsrc; \
    ann_p[NNF_Y] = ysrc;


    template<int is_descriptor, int use_allowed_patches>
    double lookup_templated(const Array<in_type> &a0, Array<double> &ann, Array<double> *ann_prev=NULL, Array<float> *coherence_temporal_sv=NULL) {
#if TABLE_VERBOSE
        printf("lookup_templated\n"); fflush(stdout);
#endif
#if TABLE_DEBUG
        check_kcoherence();
#endif

        double T0_overall = wall_time();
        get_padded_descriptor(a0, lookup_buffer, lookup_wh, "  1-NN lookup");
        Array<real> &a_wh(lookup_wh);
        
        const Array<in_type> &a(a0);
        const Array<in_type> &b(b0);
        
//        int a_w = a0.width(), a_h = a0.height();
//        int b_w = b0.width(), b_h = b0.height();

        double prop_scale = 1.0/(1+p->coherence_spatial);
        double temporal_scale = 1.0/(1+p->coherence_temporal);

        int ann_h = TABLE_IMG_HEIGHT(a_wh) - p->patch_w+1;
        int ann_w = TABLE_IMG_WIDTH(a_wh) - p->patch_w+1;
//        int b0_ew = b0.width()-p->patch_w+1;
//        int b0_eh = b0.height()-p->patch_w+1;
        int b0_ew = TABLE_IMG_WIDTH(wh0)-p->patch_w+1;
        int b0_eh = TABLE_IMG_HEIGHT(wh0)-p->patch_w+1;
    
        int a_w = ann_w+p->patch_w-1, a_h = ann_h+p->patch_w-1;
        int b_w = b0_ew+p->patch_w-1, b_h = b0_eh+p->patch_w-1;

        ann.resize(ann_h, ann_w, 3);

        PATCH_DIST_CALC_BVEC0(wh0);
        
#if TABLE_SAVE_KCOHERENCE
        Array<float> array_kcoherence_improved;
        if (p->save_kcoherence) {
            array_kcoherence_improved.resize(ann_h, ann_w, 3);
            array_kcoherence_improved.clear();
        }
#endif
        
#if TABLE_PROFILE
        double T_lookup = 0.0;
        double T_kcoherence = 0.0;
        double T_rs = 0.0;
        double T_prop = 0.0;
        
        double T0_init_inf = wall_time();
#endif
#if TABLE_VERBOSE
        printf("lookup_templated ann_prev\n"); fflush(stdout);
#endif
        if (p->init_random) {
            for (int y = 0; y < ann_h; y++) {
                double *ann_row = &ann(y, 0, 0);
                for (int x = 0; x < ann_w; x++) {
                    double *ann_p = &ann_row[x*3];
                    
                    TABLE_PRECALC_A_PATCH();
                    int xsrc = rand()%b0_ew;
                    int ysrc = rand()%b0_eh;
                    real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                    
                    ann_p[NNF_X] = xsrc;
                    ann_p[NNF_Y] = ysrc;
                    ann_p[NNF_DIST] = dcurrent;
                }
            }
        } else if (!ann_prev) {
            for (int y = 0; y < ann_h; y++) {
                double *ann_row = &ann(y, 0, 0);
                for (int x = 0; x < ann_w; x++) {
                    double *ann_p = &ann_row[x*3];
                    ann_p[NNF_X] = b0_ew+1;
                    ann_p[NNF_DIST] = DIST_INFINITY;
                }
            }
        } else {
            if (p->recalc_dist_temporal) {
                if (coherence_temporal_sv) {
                    if (coherence_temporal_sv->width() != ann_w || coherence_temporal_sv->height() != ann_h) {
                        fprintf(stderr, "coherence_temporal_sv size: %dx%d, ann size: %dx%d, a_wh size: %dx%d, a0 size: %dx%d\n", coherence_temporal_sv->width(), coherence_temporal_sv->height(), ann_w, ann_h, a_wh.width(), a_wh.height(), a0.width(), a0.height());
                        ASSERT2(false, "expected coherence_temporal_sv size to match ann size");
                    }
                }
                lookup_dist_orig.resize(ann_h, ann_w);
                for (int y = 0; y < ann_h; y++) {
                    double *ann_row = &ann(y, 0, 0);
                    double *ann_prev_row = &(*ann_prev)(y, 0, 0);
                    for (int x = 0; x < ann_w; x++) {
                        TABLE_COPY_PREV_NNF();
                        
                        TABLE_PRECALC_A_PATCH();
                        real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                        
                        double temporal_scale_current = temporal_scale;
                        if (coherence_temporal_sv) {
                            temporal_scale_current = (*coherence_temporal_sv)(y, x);
                        }
                        
                        ann_p[NNF_DIST] = lookup_dist_orig(y, x) = dcurrent * temporal_scale_current;
                    }
                }
            } else {
                for (int y = 0; y < ann_h; y++) {
                    double *ann_row = &ann(y, 0, 0);
                    double *ann_prev_row = &(*ann_prev)(y, 0, 0);
                    for (int x = 0; x < ann_w; x++) {
                        TABLE_COPY_PREV_NNF();
                        
                        ann_p[NNF_DIST] = ann_prev_p[NNF_DIST] * temporal_scale;
                    }
                }
            }
        }
#if TABLE_PROFILE
        double T_init_inf = wall_time()-T0_init_inf;
#endif
        const int dx_table[8] = {-1,  0,  1,
                                 -1,      1,
                                 -1,  0,  1};
        const int dy_table[8] = {-1, -1, -1,
                                  0,      0,
                                  1,  1,  1 };
        
        bool do_rs_or_spatial = p->do_rs || p->spatial;
        
        int prop_w = (p->query_step+p->prop_dist)*2+1;
        int sat_w = TABLE_PATCH_W+prop_w+1-1;
        int sat_w_nopad = sat_w-1;

        const int pw1 = TABLE_PATCH_W-1;
        real *wh_offset = &a_wh(pw1, pw1, 0);

#if TABLE_VERBOSE
        printf("lookup_templated begin iters\n"); fflush(stdout);
#endif
        for (int prop_iter = 0; prop_iter < p->prop_iters; prop_iter++) {
            bool reverse = prop_iter % 2;
            bool do_lookup = prop_iter == 0 && p->do_table_lookup;
            
            bool do_kcoherence = p->kcoherence > 0 && (p->kcoherence_iter < 0 || p->kcoherence_iter == prop_iter);
#if TABLE_VERBOSE
            printf("lookup_templated iter %d/%d, reverse=%d, do_kcoherence=%d\n", prop_iter, p->prop_iters, int(reverse), int(do_kcoherence)); fflush(stdout);
#endif
            
#if TABLE_OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
#endif
            for (int y0 = 0; y0 < ann_h; y0 += p->query_step) {
                int y = reverse ? (ann_h-1-y0): y0;
                
                vector<real> sat(sat_w*sat_w, 0);           /* Sat table is indexed by sat[(dy+1)*sat_w+dx+1] so we can store zeros in early rows/cols. */

                real *wh_row = wh_offset + y * a_wh.stride[0];
                for (int x0 = 0; x0 < ann_w; x0 += p->query_step) {
                    int x = reverse ? (ann_w-1-x0): x0;

#if TABLE_VERBOSE
                    printf("lookup_templated patch %d, %d, ann_w=%d, ann_h=%d\n", x, y, ann_w, ann_h); fflush(stdout);
#endif

                    TABLE_PRECALC_A_PATCH();
                    
                    int xbest = ann(y, x, NNF_X);
                    int ybest = ann(y, x, NNF_Y);
                    real dbest = ann(y, x, NNF_DIST) * prop_scale;
                    
                    /* Table lookup */

                    if (do_lookup) {
#if TABLE_VERBOSE
                        printf("lookup_templated, table lookup\n"); fflush(stdout);
#endif
#if TABLE_PROFILE
                        double T0_lookup = wall_time();
#endif
                        real *wh_p = wh_row + x * a_wh.stride[1];
#if TABLE_PRODUCT_QUANTIZE
                        if (p->product_quantize) {
                            int patch_idx = product_quantizer->quantize(wh_p);
//                            ASSERT2(in_bounds(patch_idx, table.nelems), "expected quantized patch in table bounds");
                            
                            int v = table.data[patch_idx];

                            int xsrc = INT_TO_X(v)-pw1;
                            int ysrc = INT_TO_Y(v)-pw1;

                            if (xsrc != xbest || ysrc != ybest) {
                                real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                                if (dcurrent < dbest) {
                                    dbest = dcurrent;
                                    xbest = xsrc;
                                    ybest = ysrc;
                                    CHECK_NNF_BOUNDS(x, y, xbest, ybest);
                                }
                            }
                        } else {
#endif
                            int patch_idx = 0;
                            
                            for (int i = 0; i < GRID_NCHANNELS; i++) {      /* TODO: Could speed up by templating over grid channels */
                                patch_idx += part->L[i]->get_bin(wh_p[i]) * table.stride[i];
                            }

                            int v;
#if TABLE_CLUSTER_KMEANS
                            if (p->cluster_kmeans) {
                                int v_sub = table.data[patch_idx] - sub_upper_left;
                                ASSERT2(in_bounds(v_sub, sub_tables.size()), "expected v_sub in bounds sub_tables.size()");
                                
                                PatchTableFixedN<real, in_type, ndims, itype> *sub_table = sub_tables[v_sub].get();
                                PatchPartition<real, itype> *sub_part = sub_table->part;
                                
                                int sub_patch_idx = 0;
                                for (int i = 0; i < GRID_NCHANNELS; i++) {
                                    sub_patch_idx += sub_part->L[i]->get_bin(wh_p[i]) * sub_table->table.stride[i];
                                }
                                ASSERT2(in_bounds(sub_patch_idx, sub_table->table.nelems), "expected sub_patch_idx in bounds for sub_table");
                                v = sub_table->table.data[sub_patch_idx];
                            } else {
                                v = table.data[patch_idx];
                            }
#else
                            v = table.data[patch_idx];
#endif

                            int xsrc = INT_TO_X(v)-pw1;
                            int ysrc = INT_TO_Y(v)-pw1;

                            if (xsrc != xbest || ysrc != ybest) {
                                real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                                if (dcurrent < dbest) {
                                    dbest = dcurrent;
                                    xbest = xsrc;
                                    ybest = ysrc;
                                    CHECK_NNF_BOUNDS(x, y, xbest, ybest);
                                }
                            }
                            
#if TABLE_MULTI_TABLES
                            int xprev = xsrc, yprev = ysrc;
                            for (int itable = 0; itable < (int) tables.size(); itable++) {
                                auto current_table = tables[itable].get();
                                
                                patch_idx = 0;
                                for (int i = 0; i < GRID_NCHANNELS; i++) {      /* TODO: Could speed up by templating over grid channels */
                                    patch_idx += current_table->part->L[i]->get_bin(wh_p[i]) * current_table->table.stride[i];
                                }
                                
                                v = current_table->table.data[patch_idx];
                                
                                int xsrc1 = INT_TO_X(v)-pw1;
                                int ysrc1 = INT_TO_Y(v)-pw1;
                                
                                if ((xsrc1 != xbest || ysrc1 != ybest) && (xsrc1 != xprev || ysrc1 != yprev)) {
                                    real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc1, ysrc1);
                                    if (dcurrent < dbest) {
                                        dbest = dcurrent;
                                        xbest = xsrc1;
                                        ybest = ysrc1;
                                        CHECK_NNF_BOUNDS(x, y, xbest, ybest);
                                    }
                                }
                                xprev = xsrc1;
                                yprev = ysrc1;
                            }
#endif
#if TABLE_PRODUCT_QUANTIZE
                        }
#endif
                        
#if TABLE_PROFILE
                        T_lookup += wall_time() - T0_lookup;
#endif
                    }
                    
                    /* k-coherence */
                    
                    if (do_kcoherence) {
#if TABLE_VERBOSE
                        printf("lookup_templated, k-coherence, xbest=%d, ybest=%d, kcoherence_set, width=%d, height=%d, kcoherence_step=%d\n", xbest, ybest, kcoherence_set.width(), kcoherence_set.height(), p->kcoherence_step); fflush(stdout);
#endif
#if TABLE_PROFILE
                        double T0_kcoherence = wall_time();
#endif
                        /* TODO: Does it work better to use the kcoherence set in the new position if (xbest, ybest) changes? */
#if TABLE_KCOHERENCE_STEP
                        kcoherence_type *kcoherence_ptr = kcoherence_set.data + kcoherence_set.stride[0]*(ybest/p->kcoherence_step) + kcoherence_set.stride[1]*(xbest/p->kcoherence_step);
#if TABLE_DEBUG
                        if (allowed_patches && (*allowed_patches)(ybest/p->kcoherence_step, xbest/p->kcoherence_step) != p->allowed_index) {
                            fprintf(stderr, "kcoherence accessed disallowed patch %d, %d (xbest=%d, ybest=%d, kcoherence_step=%d)\n", xbest/p->kcoherence_step, ybest/p->kcoherence_step, p->kcoherence_step); ASSERT2(false, "kcoherence access disallowed patch");
                        }
#endif
                        
                        int kcoherence_dx = 0, kcoherence_dy = 0;
                        if (p->kcoherence_step > 1) {
                            kcoherence_dx = xbest%p->kcoherence_step;
                            kcoherence_dy = ybest%p->kcoherence_step;
                        }
#else
                        kcoherence_type *kcoherence_ptr = kcoherence_set.data + kcoherence_set.stride[0]*(ybest) + kcoherence_set.stride[1]*(xbest);
#endif
                        for (int k = 0; k < p->kcoherence; k++) {
#if TABLE_VERBOSE
                            printf("lookup_templated, kcoherence %d/%d\n", k, p->kcoherence);
#endif
                            kcoherence_type *k_p = &kcoherence_ptr[k];
#if TABLE_KCOHERENCE_TRIANGLE
                            float k_dist = TABLE_KCOHERENCE_DIST(*k_p);
                            if (k_dist >= p->triangle_factor*dbest) {
                                break;
                            }
#endif
                            int v = TABLE_KCOHERENCE_POS(*k_p);
                            int xsrc = INT_TO_X(v);
                            int ysrc = INT_TO_Y(v);
#if TABLE_KCOHERENCE_STEP
                            if (p->kcoherence_step > 1) {
                                xsrc += kcoherence_dx;
                                ysrc += kcoherence_dy;
                                if (xsrc >= b0_ew) { xsrc = b0_ew-1; }
                                if (ysrc >= b0_eh) { ysrc = b0_eh-1; }
                            }
#endif
                            CHECK_NNF_BOUNDS(x, y, xsrc, ysrc);
                            real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);

                            if (dcurrent < dbest) {
#if TABLE_SAVE_KCOHERENCE
                                if (p->save_kcoherence) {
                                    for (int k_channel = 0; k_channel < 3; k_channel++) {
                                        array_kcoherence_improved.get_nearest(y, x, k_channel) += 1.0/(p->kcoherence*p->prop_iters);
                                    }
                                }
#endif

                                
                                dbest = dcurrent;
                                xbest = xsrc;
                                ybest = ysrc;
                            }
                            
#if TABLE_KCOHERENCE_ENRICH
                            if (p->kcoherence_enrich) {
                                kcoherence_type *kcoherence_ptr2 = kcoherence_set.data + kcoherence_set.stride[0]*(ysrc/p->kcoherence_step) + kcoherence_set.stride[1]*(xsrc/p->kcoherence_step);
                                int kcoherence_dx2 = 0, kcoherence_dy2 = 0;
                                if (p->kcoherence_step > 1) {
                                    kcoherence_dx2 = xsrc%p->kcoherence_step;
                                    kcoherence_dy2 = ysrc%p->kcoherence_step;
                                }
                                for (int k2 = 0; k2 < p->kcoherence; k2++) {
                                    int v2 = TABLE_KCOHERENCE_POS(kcoherence_ptr2[k]);
                                    int xsrc2 = INT_TO_X(v2);
                                    int ysrc2 = INT_TO_Y(v2);
                                    if (p->kcoherence_step > 1) {
                                        xsrc2 += kcoherence_dx2;
                                        ysrc2 += kcoherence_dy2;
                                        if (xsrc2 >= b0_ew) { xsrc2 = b0_ew-1; }
                                        if (ysrc2 >= b0_eh) { ysrc2 = b0_eh-1; }
                                    }
                                    CHECK_NNF_BOUNDS(x, y, xsrc2, ysrc2);

                                    if (xsrc2 != xbest || ysrc2 != ybest) {
                                        real dcurrent2 = TABLE_LOOKUP_PATCH_DIST(xsrc2, ysrc2);
                                        
                                        if (dcurrent2 < dbest) {
                                            dbest = dcurrent2;
                                            xbest = xsrc2;
                                            ybest = ysrc2;
                                        }
                                    }

                                }
                            }
#endif
                        }
#if TABLE_PROFILE
                        T_kcoherence += wall_time() - T0_kcoherence;
#endif
                    }
                    
                    /* Spatial search or random search */
                    if (do_rs_or_spatial) {
#if TABLE_VERBOSE
                        printf("lookup_templated, rs\n"); fflush(stdout);
#endif
#if TABLE_PROFILE
                        double T0_rs = wall_time();
#endif
                        if (p->spatial) {
                            int spatial_ymin = ybest-p->spatial;
                            int spatial_ymax = ybest+p->spatial+1;
                            int spatial_xmin = xbest-p->spatial;
                            int spatial_xmax = xbest+p->spatial+1;
                            int xbest0 = xbest, ybest0 = ybest;
                            if (spatial_ymin < 0) { spatial_ymin = 0; }
                            else if (spatial_ymax > b0_eh) { spatial_ymax = b0_eh; }
                            if (spatial_xmin < 0) { spatial_xmin = 0; }
                            else if (spatial_xmax > b0_ew) { spatial_xmax = b0_ew; }
                            
                            for (int ysrc = spatial_ymin; ysrc < spatial_ymax; ysrc++) {
                                for (int xsrc = spatial_xmin; xsrc < spatial_xmax; xsrc++) {
                                    if ((xsrc != xbest0 || ysrc != ybest0) && (!use_allowed_patches || (*allowed_patches)(ysrc, xsrc) == p->allowed_index)) {
                                        real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                                        
                                        if (dcurrent < dbest) {
                                            dbest = dcurrent;
                                            xbest = xsrc;
                                            ybest = ysrc;
                                            CHECK_NNF_BOUNDS(x, y, xsrc, ysrc);
                                        }
                                    }
                                }
                            }
                        } else {
                            int r = rand()&7;
                            int xsrc = xbest+dx_table[r];
                            int ysrc = ybest+dy_table[r];
                            
                            if (in_bounds(xsrc, b0_ew) && in_bounds(ysrc, b0_eh) && (!use_allowed_patches || (*allowed_patches)(ysrc, xsrc) == p->allowed_index)) {
                                real dcurrent = TABLE_LOOKUP_PATCH_DIST(xsrc, ysrc);
                                
                                if (dcurrent < dbest) {
                                    dbest = dcurrent;
                                    xbest = xsrc;
                                    ybest = ysrc;
                                    CHECK_NNF_BOUNDS(x, y, xsrc, ysrc);
                                }
                            }
                        }
#if TABLE_PROFILE
                        T_rs += wall_time() - T0_rs;
#endif
                    }
                
                    /* Propagate */
                    if (p->do_prop) {
#if TABLE_VERBOSE
                        printf("lookup_templated, propagate\n"); fflush(stdout);
#endif
#if TABLE_PROFILE
                        double T0_prop = wall_time();
#endif

                        int ax_shift = x - p->query_step;
                        int ay_shift = y - p->query_step;
                        int bx_shift = xbest - p->query_step;
                        int by_shift = ybest - p->query_step;
                        real *a_base = a.data + ay_shift*a.stride[0] + ax_shift*a.stride[1];
                        real *b_base = b.data + by_shift*b.stride[0] + bx_shift*b.stride[1];
                        double *ann_base = ann.data + ay_shift*ann.stride[0] + ax_shift*ann.stride[1];
                        
                        int dy_min = 0, dy_max = sat_w_nopad;
                        int dx_min = 0, dx_max = sat_w_nopad;

                        if (ax_shift + dx_min < 0) { dx_min = -ax_shift; }
                        else if (ax_shift + dx_max > a_w) { dx_max = a_w-ax_shift; }
                        if (ay_shift + dy_min < 0) { dy_min = -ay_shift; }
                        else if (ay_shift + dy_max > a_h) { dy_max = a_h-ay_shift; }
                        
                        if (bx_shift + dx_min < 0) { dx_min = -bx_shift; }
                        else if (bx_shift + dx_max > b_w) { dx_max = b_w-bx_shift; }
                        if (by_shift + dy_min < 0) { dy_min = -by_shift; }
                        else if (by_shift + dy_max > b_h) { dy_max = b_h-by_shift; }

                        ASSERT(dx_min >= 0, "expected dx_min >= 0");
                        ASSERT(dy_min >= 0, "expected dy_min >= 0");
                        ASSERT(dx_max <= sat_w_nopad, "expected dx_max <= sat_w_nopad");
                        ASSERT(dy_max <= sat_w_nopad, "expected dy_max <= sat_w_nopad");

#if TABLE_PROP_FAST
                        if (!is_descriptor) {
#if TABLE_VERBOSE
                            printf("lookup_templated, not is_descriptor, propagating, zero sat\n"); fflush(stdout);
#endif
                            if (dx_min > 0 || dy_min > 0) {
                                if (dx_min > 0) {
                                    /* Zero column left of dx_min in sat table */
                                    for (int dy = dy_min > 0 ? dy_min-1: dy_min; dy < dy_max; dy++) {
                                        sat[(dy+1)*sat_w+dx_min] = 0;
                                    }
                                }
                                if (dy_min > 0) {
                                    /* Zero row above dy_min in sat table */
                                    for (int dx = dx_min; dx < dx_max; dx++) {
                                        sat[(dy_min)*sat_w+dx+1] = 0;
                                    }
                                }
                            }
#if TABLE_VERBOSE
                            printf("lookup_templated, not is_descriptor, propagating, compute sat\n"); fflush(stdout);
#endif
                            
                            /* Compute summed area table */
                            for (int dy = dy_min; dy < dy_max; dy++) {
                                real *a_row = a_base + dy*a.stride[0];
                                real *b_row = b_base + dy*b.stride[0];
                                real *sat_row = &sat[(dy+1)*sat_w+1];
                                real *sat_row_prev = &sat[dy*sat_w+1];
                                real row_sum = 0;
                                
                                for (int dx = dx_min; dx < dx_max; dx++) {            /* TODO: Could probably unroll this via a template for speed */
                                    ASSERT(in_bounds(ax_shift+dx, a_w), "expected ax_shift+dx in bounds a_w");
                                    ASSERT(in_bounds(ay_shift+dy, a_h), "expected ay_shift+dy in bounds a_h");
                                    ASSERT(in_bounds(bx_shift+dx, b_w), "expected bx_shift+dx in bounds b_w");
                                    ASSERT(in_bounds(by_shift+dy, b_h), "expected by_shift+dy in bounds b_h");
                                    
                                    real *a_pixel = a_row + dx*a.stride[1];
                                    real *b_pixel = b_row + dx*b.stride[1];
                                    
                                    real delta_0 = a_pixel[0]-b_pixel[0];
                                    real delta_1 = a_pixel[1]-b_pixel[1];
                                    real delta_2 = a_pixel[2]-b_pixel[2];
                                    
                                    row_sum += delta_0*delta_0 + delta_1*delta_1 + delta_2*delta_2;
                                    sat_row[dx] = row_sum + sat_row_prev[dx];
                                }
                            }
                        }
#endif
                        
#if TABLE_VERBOSE
                        printf("lookup_templated, propagating, result from sat\n"); fflush(stdout);
#endif

                        /* Compute patch distances. TODO: No need to store for central patch. Also can we compute it faster knowing central patch dist? */
                        int dy_max_patch = dy_max - TABLE_PATCH_W+1;
                        int dx_max_patch = dx_max - TABLE_PATCH_W+1;
    //                    if (x == 100 && y == 100) { printf("dy: %d %d, dx: %d %d\n", dy_min, dy_max_patch, dx_min, dx_max_patch); }
    /*                    if (x == 1324 && y == 0) {
                            printf("at %d, %d, dx_min=%d, dx_max=%d, dy_min=%d, dy_max=%d, xbest=%d, ybest=%d, dbest=%f, dx_max_patch=%d, dy_max_patch=%d\n", x, y, dx_min, dx_max, dy_min, dy_max, xbest, ybest, dbest, dx_max_patch, dy_max_patch);
                        } */
                        
                        bool missing = false;
                        for (int dy = dy_min; dy < dy_max_patch; dy++) {
    //                    for (int dy = dy_max_patch - 1; dy >= dy_min; dy--) {
#if TABLE_PROP_FAST
                            real *sat_row_1;
                            real *sat_row_2;
                            if (!is_descriptor) {
                                sat_row_1 = &sat[(dy)*sat_w];
                                sat_row_2 = &sat[(dy+TABLE_PATCH_W)*sat_w];
                            }
#endif
                            double *ann_row = ann_base + dy*ann.stride[0];
                            int ysrc = dy + by_shift;
                            itype *allowed_patches_row;
                            if (use_allowed_patches) {
                                allowed_patches_row = allowed_patches->data + ysrc * allowed_patches->stride[0];
                            }
                            ASSERT(in_bounds(ysrc, b0_eh), "expected ysrc in bounds in propagate");

                            for (int dx = dx_min; dx < dx_max_patch; dx++) {
    //                        for (int dx = dx_max_patch-1; dx >= dx_min; dx--) {
#if TABLE_VERBOSE
                                if (!in_bounds(ax_shift+dx, ann_w) || !in_bounds(ay_shift+dy, ann_h)) {
                                    fprintf(stderr, "x=%d, y=%d, dx=%d, dy=%d, ax_shift=%d, ay_shift=%d, ann_w=%d, ann_h=%d, a_w=%d, a_h=%d, dx_min=%d, dx_max=%d, dx_max_patch=%d, dy_min=%d, dy_max=%d, dy_max_patch=%d, bx_shift=%d, by_shift=%d, query_step=%d, sat_w_nopad=%d, b_w=%d, b_h=%d, negative one: %d\n",
                                            x, y, dx, dy, ax_shift, ay_shift, ann_w, ann_h, a_w, a_h, dx_min, dx_max, dx_max_patch, dy_min, dy_max, dy_max_patch, bx_shift, by_shift, p->query_step, sat_w_nopad, b_w, b_h, -1); ASSERT2(false, "out of bounds");
                                }
#endif
                                ASSERT(in_bounds(ax_shift+dx+TABLE_PATCH_W-1, a_w), "expected ax_shift+dx+patch_w-1 in bounds a_w");
                                ASSERT(in_bounds(ay_shift+dy+TABLE_PATCH_W-1, a_h), "expected ay_shift+dy+patch_w-1 in bounds a_h");
                                ASSERT(in_bounds(ax_shift+dx, ann_w), "expected ax_shift+dx in bounds ann_w");
                                ASSERT(in_bounds(ay_shift+dy, ann_h), "expected ay_shift+dy in bounds ann_h");
                                real dcurrent;
                                int xsrc = dx + bx_shift;
                                ASSERT(in_bounds(xsrc, b0_ew), "expected xsrc in bounds in propagate");
                                if (use_allowed_patches && allowed_patches_row[xsrc] != p->allowed_index) {
                                    missing = true;
                                    continue;
                                }
                                if (!is_descriptor) {
#if TABLE_PROP_FAST
                                    int dx2 = dx + TABLE_PATCH_W;
                                    dcurrent = sat_row_1[dx] - sat_row_1[dx2] - sat_row_2[dx] + sat_row_2[dx2];
#else
                                    real *avec_shifted = &a.get_nearest(dy+ay_shift, dx+ax_shift, 0);
                                    dcurrent = patch_dist_exact<real>(a, avec_shifted, b0, xsrc, ysrc);
#endif
                                } else {
                                    real *avec_shifted = &a_wh.get_nearest(dy+ay_shift+pw1, dx+ax_shift+pw1, 0);
                                    dcurrent = patch_dist_approx<real, TABLE_NCHANNELS>(wh0, avec_shifted, bvec0, xsrc, ysrc);
                                }

                                double *ann_ptr = ann_row + dx*ann.stride[1];
                                if (dcurrent < ann_ptr[NNF_DIST]) {
#if TABLE_DEBUG
                                    if (!is_descriptor) {
                                        real *avec_shifted = &a.get_nearest(dy+ay_shift, dx+ax_shift, 0);
                                        real d_recalc = patch_dist_exact<real>(a, avec_shifted, b0, dx+bx_shift, dy+by_shift);
                                        if (fabs(d_recalc-dcurrent) > 1e-3) {
                                            fprintf(stderr, "dcurrent=%f, d_recalc=%f, differ, dx=%d, dy=%d, x, y=%d, %d, xbest=%d, ybest=%d, target b position: %d, %d, dx_min, dx_max: %d %d, dy_min, dy_max: %d %d, sat_w_nopad: %d, ax_shift: %d, ay_shift: %d, bx_shift: %d, by_shift: %d, dx_max_patch: %d, dy_max_patch: %d, negative 1: %d\n", dcurrent, d_recalc, dx, dy, x, y, xbest, ybest, dx+bx_shift, dy+by_shift, dx_min, dx_max, dy_min, dy_max, sat_w_nopad, ax_shift, ay_shift, bx_shift, by_shift, dx_max_patch, dy_max_patch, -1);
                                            printf("sat table:\n");
                                            for (int iy = 0; iy < sat_w; iy++) {
                                                for (int ix = 0; ix < sat_w; ix++) {
                                                    printf("%.4f ", sat[iy*sat_w+ix]);
                                                }
                                                printf("\n");
                                            }
                                            printf("\n");
                                            ASSERT2(false, "differ");
                                        }
                                    }
                                    //fprintf(stderr, "success: dcurrent=%f, d_recalc=%f, differ, dx=%d, dy=%d, x, y=%d, %d, xbest=%d, ybest=%d, target b position: %d, %d, dx_min, dx_max: %d %d, dy_min, dy_max: %d %d, sat_w_nopad: %d, ax_shift: %d, ay_shift: %d, bx_shift: %d, by_shift: %d, dx_max_patch: %d, dy_max_patch: %d, negative 1: %d\n", dcurrent, d_recalc, dx, dy, x, y, xbest, ybest, dx+bx_shift, dy+by_shift, dx_min, dx_max, dy_min, dy_max, sat_w_nopad, ax_shift, ay_shift, bx_shift, by_shift, dx_max_patch, dy_max_patch, -1);
#endif
                                    ann_ptr[NNF_X] = xsrc;
                                    ann_ptr[NNF_Y] = ysrc;
                                    ann_ptr[NNF_DIST] = dcurrent;
                                }
                            }
                        }
#if TABLE_VERBOSE
                        printf("lookup_templated, propagating, filling missing\n"); fflush(stdout);
#endif
                        if (dx_max < sat_w_nopad || dy_max < sat_w_nopad || missing) {
                            int dy_start, dy_end;
                            int dx_start, dx_end;
                            if (missing) {
                                dx_start = 0; dx_end = p->query_step*2;
                                dy_start = 0; dy_end = p->query_step*2;
                                if (ax_shift + dx_start < 0) { dx_start = -ax_shift; }
                                if (ay_shift + dy_start < 0) { dy_start = -ay_shift; }
                                if (ax_shift + dx_end > ann_w) { dx_end = ann_w-ax_shift; }
                                if (ay_shift + dy_end > ann_h) { dy_end = ann_h-ay_shift; }
                            } else {
                                if (reverse) {
                                    dx_start = 0; dx_end = p->query_step;
                                    dy_start = 0; dy_end = p->query_step;
                                    if (ax_shift + dx_start < 0) { dx_start = -ax_shift; }
                                    if (ay_shift + dy_start < 0) { dy_start = -ay_shift; }
                                } else {
                                    dx_start = p->query_step; dx_end = p->query_step*2;
                                    dy_start = p->query_step; dy_end = p->query_step*2;
                                    if (ax_shift + dx_end > ann_w) { dx_end = ann_w-ax_shift; }
                                    if (ay_shift + dy_end > ann_h) { dy_end = ann_h-ay_shift; }
                                }
                            }

                            for (int dy = dy_start; dy < dy_end; dy++) {
                                double *ann_row = ann_base + dy*ann.stride[0];
                                int ysrc = ybest + (ay_shift+dy)-(y);
                                if (ysrc < 0) { ysrc = 0; }
                                else if (ysrc >= b0_eh) { ysrc = b0_eh-1; }

                                for (int dx = dx_start; dx < dx_end; dx++) {
                                    double *ann_ptr = ann_row + dx*ann.stride[1];
                                    if (ann_ptr[NNF_DIST] == DIST_INFINITY) {
                                        int xsrc = xbest + (ax_shift+dx)-(x);
                                        if (xsrc < 0) { xsrc = 0; }
                                        else if (xsrc >= b0_ew) { xsrc = b0_ew-1; }

                                        if (use_allowed_patches) {
                                            while ((*allowed_patches)(ysrc, xsrc) != p->allowed_index) {
                                                xsrc = rand()%b0_ew;
                                                ysrc = rand()%b0_eh;
                                            }
                                        }

                                        real dcurrent;
                                        if (is_descriptor) {
                                            real *avec_shifted = &a_wh.get_nearest(dy+ay_shift+pw1, dx+ax_shift+pw1, 0);
                                            dcurrent = patch_dist_approx<real, TABLE_NCHANNELS>(wh0, avec_shifted, bvec0, xsrc, ysrc);
                                        } else {
                                            real *avec_shifted = &a.get_nearest(dy+ay_shift, dx+ax_shift, 0);
                                            dcurrent = patch_dist_exact<real>(a, avec_shifted, b0, xsrc, ysrc);
                                        }
                                        ann_ptr[NNF_DIST] = dcurrent;
                                        ann_ptr[NNF_X] = xsrc;
                                        ann_ptr[NNF_Y] = ysrc;
                                    }
                                }
                            }

                        }
#if TABLE_PROFILE
                        T_prop += wall_time() - T0_prop;
#endif
                    } else {
                        ann(y, x, NNF_X) = xbest;
                        ann(y, x, NNF_Y) = ybest;
                        ann(y, x, NNF_DIST) = dbest;
                    }
//#if TABLE_DEBUG
//                    fprintf(stderr, "end successful loop\n\n");
//#endif
                    //if (dx_min > 0 || dy_min > 0 || dx_max < sat_w_nopad || dy_max < sat_w_nopad) {
                    //    if (dx_min > 0) {
                    //
                    //    }
                    //}
                    /*
                    ann(y, x, NNF_X) = xbest;
                    ann(y, x, NNF_Y) = ybest;
                    ann(y, x, NNF_DIST) = dbest;
                    if (fabs(ann(y, x, NNF_X) - xbest) > 1e-5 || fabs(ann(y, x, NNF_Y) - ybest) > 1e-5 || fabs(ann(y, x, NNF_DIST)-dbest) > 1e-5) {
                        fprintf(stderr, "expected ann(y, x, :) to match best: %f %f %f, %f %f %f\n", ann(y, x, NNF_X), ann(y, x, NNF_Y), ann(y, x, NNF_DIST),
                        double(xbest), double(ybest), dbest);
                        ASSERT2(false, "failed");
                    }
                    */
                }
            }
            
            
        }

        if (ann_prev && p->coherence_temporal != 0.0) {
            for (int y = 0; y < ann_h; y++) {
                double *ann_row = &ann(y, 0, 0);
                double *ann_prev_row = &(*ann_prev)(y, 0, 0);
                for (int x = 0; x < ann_w; x++) {
                    double *ann_p = &ann_row[x*3];
                    double *ann_prev_p = &ann_prev_row[x*3];
                    if (ann_p[NNF_X] == ann_prev_p[NNF_X] && ann_p[NNF_Y] == ann_prev_p[NNF_Y]) {
                        if ((p->recalc_dist_temporal && ann_p[NNF_DIST] == lookup_dist_orig(y, x)) ||
                            (!p->recalc_dist_temporal && ann_p[NNF_DIST] == ann_prev_p[NNF_DIST] * temporal_scale)) {
                            ann_p[NNF_DIST] /= temporal_scale;
                        }
                    }
                }
            }
        }
        
#if TABLE_SAVE_KCOHERENCE
        if (p->save_kcoherence) {
            save_color_image<float>(array_kcoherence_improved, "kcoherence_improved.pfm");
        }
#endif


#if TABLE_PROFILE
        printf("Profile times:\n"); fflush(stdout);
        double T_sum = T_init_inf + T_lookup + T_kcoherence + T_rs + T_prop;
        printf("  T_init_inf:   %7.4f sec (%5.2f%%)\n", T_init_inf, 100.0*T_init_inf/T_sum);
        printf("  T_lookup:     %7.4f sec (%5.2f%%)\n", T_lookup, 100.0*T_lookup/T_sum);
        printf("  T_kcoherence: %7.4f sec (%5.2f%%)\n", T_kcoherence, 100.0*T_kcoherence/T_sum);
        printf("  T_rs:         %7.4f sec (%5.2f%%)\n", T_rs, 100.0*T_rs/T_sum);
        printf("  T_prop:       %7.4f sec (%5.2f%%)\n", T_prop, 100.0*T_prop/T_sum);
#endif
        double ans = wall_time()-T0_overall;
        if (p->verbose) {
            printf("    1-NN overall lookup time: %f secs\n", ans);
        }

        return ans;
    }
    
    double lookup(const Array<in_type> &a0, Array<double> &ann, Array<double> *ann_prev=NULL, Array<float> *coherence_temporal_sv=NULL) {
#if TABLE_VERBOSE
        printf("lookup\n"); fflush(stdout);
#endif
        if (p->lookup_algo == LOOKUP_ALGO_TREECANN) {
            double T0 = wall_time();
            get_padded_descriptor(a0, lookup_buffer, lookup_wh, "  1-NN lookup");
            treecann->lookup(a0, lookup_wh, ann);
            return wall_time() - T0;
        }
        
        if (p->is_descriptor) {
            if (allowed_patches) {
                return lookup_templated<1, 1>(a0, ann, ann_prev, coherence_temporal_sv);
            } else {
                return lookup_templated<1, 0>(a0, ann, ann_prev, coherence_temporal_sv);
            }
        } else {
            if (allowed_patches) {
                return lookup_templated<0, 1>(a0, ann, ann_prev, coherence_temporal_sv);
            } else {
                return lookup_templated<0, 0>(a0, ann, ann_prev, coherence_temporal_sv);
            }
        }
    }
	void euclidean_dt(float *fuction, int n, float *distance, int *current_parabola, float* previous_intersection, int* index_array) {
		int preNo = 0;
		current_parabola[0] = 0;
		previous_intersection[0] = -INF;
		previous_intersection[1] = +INF;
		for (int curNo = 1; curNo <= n - 1; curNo++) {
			float  cur_intersection = ((fuction[curNo] + curNo*curNo) - (fuction[current_parabola[preNo]] + current_parabola[preNo] * current_parabola[preNo])) / (2 * curNo - 2 * current_parabola[preNo]);
			while (cur_intersection <= previous_intersection[preNo]) {
				preNo--;
				cur_intersection = ((fuction[curNo] + curNo*curNo) - (fuction[current_parabola[preNo]] + current_parabola[preNo] * current_parabola[preNo])) / (2 * curNo - 2 * current_parabola[preNo]);
			}
			preNo++;
			current_parabola[preNo] = curNo;
			previous_intersection[preNo] = cur_intersection;
			previous_intersection[preNo + 1] = +INF;
		}

		preNo = 0;
		for (int curNo = 0; curNo <= n - 1; curNo++) {
			while (previous_intersection[preNo + 1] < curNo)
				preNo++;
			distance[curNo] = (curNo - current_parabola[preNo])*(curNo - current_parabola[preNo]) + fuction[current_parabola[preNo]];
			// Connelly: Store patch index in index_array[q]
			index_array[curNo] = index_array[current_parabola[preNo]];
		}
	}
	//add by liming April 2nd.2015
    // Euclidean distance transform algorithm
    // Felzenszwalb and Huttenlocher, Distance Transforms of Sampled Functions, Theory of Computing 2012
	void do_dt_euclidean() {
        double T_begin_dt = wall_time();
		Array<float> table_dist(table.sizes);
		for (int j = 0; j < table.nelems; j++)
		{
			if (table.data[j] == TABLE_UNUSED)
			{
				table_dist.data[j] = INF;
			}
			else
			{
				table_dist.data[j] = 0.0;
			}
		}
		// Connelly: Create table_dist and initialize the values in it:
		// vector<float> table_dist(table.sizes);
		// for (int j = 0; j < table.nelems; j++) {
		//   // Do an if test to initialize distances
		// }
		int maxdim = 0;
		for (int n = 0; n < table_dist.sizes.size(); n++)
		{
			maxdim = std::max(maxdim, table_dist.sizes[n]);
		}
		float *fuction = new float[maxdim];
		int *current_parabola = new int[maxdim];
		float *previous_intersection = new float[maxdim + 1];
		float *distance = new float[maxdim];
		int * index_array = new int[maxdim];
        
		vector<int> outer_stride(table_dist.sizes.size());
		for (int inner_dim = 0; inner_dim < table_dist.sizes.size(); inner_dim++) {		// Dimension to be transformed along
			int product = 1;		// Product of all the rest of the dimensions (that are not the inner dimension)
			for (int dim = table_dist.sizes.size() - 1; dim >= 0; dim--) {
				if (dim != inner_dim) {
					outer_stride[dim] = product;
					product *= table_dist.sizes[dim];
				}
			}

			for (int j = 0; j < product; j++) {		// Index within the product of all the remaining dimensions
				// Convert index j into correct index of array
				int array_index = 0;
				for (int dim = 0; dim < table_dist.sizes.size(); dim++) {
					if (dim != inner_dim) {
						// Decode index j into the index within the given dimension
						int index_within_dim = (j / outer_stride[dim]) % table_dist.sizes[dim];
						array_index += index_within_dim * table_dist.stride[dim];
					}
				}

				// (in PatchTableFixedN  itype *index_array    -- because itype is the integer type)

				int stride_inner_dim = table_dist.stride[inner_dim];
				for (int inner_index = 0; inner_index < table_dist.sizes[inner_dim]; inner_index++) {
					fuction[inner_index] = table_dist.data[array_index + inner_index * stride_inner_dim];
					// Connelly: Write to index_array by reading from table.data[...] where ... is the same thing as above;
					index_array[inner_index] = table.data[array_index + inner_index * stride_inner_dim];
				}

				// Connelly: Call dt(f, table_dist.sizes[inner_dim], d, index_array, v, z);   -- More arguments to dt
				euclidean_dt(fuction, table_dist.sizes[inner_dim], distance, current_parabola, previous_intersection, index_array);

				for (int inner_index = 0; inner_index < table_dist.sizes[inner_dim]; inner_index++) {
					table_dist.data[array_index + inner_index * stride_inner_dim] = distance[inner_index];
					// Connelly: Write back the computed index to 'table'
					table.data[array_index + inner_index * stride_inner_dim] = index_array[inner_index];
				}
			}
		}
        delete[] index_array;
		delete[] distance;
		delete[] fuction;
		delete[] current_parabola;
		delete[] previous_intersection;
        if (p->verbose) {
            double T_end_dt = wall_time();
            printf("table dt time (euclidean): %f secs\n", T_end_dt-T_begin_dt);
        }
	}
};

/* ------------------------------------------------------------------------------------------------
   PatchTable class with variable dimension count
   ------------------------------------------------------------------------------------------------ */

template<class real=float, class in_type=float, class itype=int32_t>
class PatchTable { public:
    int gck_dims;
    int gck_nchroma;

#define DECLARE_PTR(dims) \
    PatchTableFixedN<real, in_type, dims, itype> *p##dims;

#define DECLARE_CONSTRUCTOR(dims) \
        else if (gck_dims == dims) { p##dims = new PatchTableFixedN<real, in_type, dims, itype>(p, a, allowed_patches); }

#define DECLARE_DESTRUCTOR(dims) \
        else if (gck_dims == dims) { delete p##dims; }

#define DECLARE_LOOKUP(dims) \
        else if (gck_dims == dims) { return p##dims->lookup(a, ann, ann_prev, coherence_temporal_sv); }

#define DECLARE_SAVE(dims) \
        else if (gck_dims == dims) { p##dims->save(filename); }

    DECLARE_PTR(3)
    DECLARE_PTR(4)
    DECLARE_PTR(5)
    DECLARE_PTR(6)
    DECLARE_PTR(7)
    DECLARE_PTR(8)
    DECLARE_PTR(9)
    DECLARE_PTR(10)
    DECLARE_PTR(11)
    DECLARE_PTR(12)
    DECLARE_PTR(13)
    DECLARE_PTR(14)
    DECLARE_PTR(15)
    DECLARE_PTR(20)
    DECLARE_PTR(25)
    DECLARE_PTR(30)
    DECLARE_PTR(40)
    
    PatchTable(PatchTableParams *p, const Array<in_type> &a, Array<itype> *allowed_patches=NULL) {
        if (p->is_descriptor) {
            p->ndims = a.channels();
            p->nchroma = 0;
        }
        gck_nchroma = p->nchroma;
        gck_dims = p->ndims;
        
        if (0) { }
        DECLARE_CONSTRUCTOR(3)
        DECLARE_CONSTRUCTOR(4)
        DECLARE_CONSTRUCTOR(5)
        DECLARE_CONSTRUCTOR(6)
        DECLARE_CONSTRUCTOR(7)
        DECLARE_CONSTRUCTOR(8)
        DECLARE_CONSTRUCTOR(9)
        DECLARE_CONSTRUCTOR(10)
        DECLARE_CONSTRUCTOR(11)
        DECLARE_CONSTRUCTOR(12)
        DECLARE_CONSTRUCTOR(13)
        DECLARE_CONSTRUCTOR(14)
        DECLARE_CONSTRUCTOR(15)
        DECLARE_CONSTRUCTOR(20)
        DECLARE_CONSTRUCTOR(25)
        DECLARE_CONSTRUCTOR(30)
        DECLARE_CONSTRUCTOR(40)
        else { fprintf(stderr, "unsupported number of dimensions\n"); ASSERT2(false, "invalid dims"); }
    }
    
    ~PatchTable() {
        if (0) { }
        DECLARE_DESTRUCTOR(3)
        DECLARE_DESTRUCTOR(4)
        DECLARE_DESTRUCTOR(5)
        DECLARE_DESTRUCTOR(6)
        DECLARE_DESTRUCTOR(7)
        DECLARE_DESTRUCTOR(8)
        DECLARE_DESTRUCTOR(9)
        DECLARE_DESTRUCTOR(10)
        DECLARE_DESTRUCTOR(11)
        DECLARE_DESTRUCTOR(12)
        DECLARE_DESTRUCTOR(13)
        DECLARE_DESTRUCTOR(14)
        DECLARE_DESTRUCTOR(15)
        DECLARE_DESTRUCTOR(20)
        DECLARE_DESTRUCTOR(25)
        DECLARE_DESTRUCTOR(30)
        DECLARE_DESTRUCTOR(40)
        else { fprintf(stderr, "unsupported number of dimensions\n"); ASSERT2(false, "invalid dims"); }
    }

    double lookup(const Array<in_type> &a, Array<double> &ann, Array<double> *ann_prev=NULL, Array<float> *coherence_temporal_sv=NULL) {
        if (0) { }
        DECLARE_LOOKUP(3)
        DECLARE_LOOKUP(4)
        DECLARE_LOOKUP(5)
        DECLARE_LOOKUP(6)
        DECLARE_LOOKUP(7)
        DECLARE_LOOKUP(8)
        DECLARE_LOOKUP(9)
        DECLARE_LOOKUP(10)
        DECLARE_LOOKUP(11)
        DECLARE_LOOKUP(12)
        DECLARE_LOOKUP(13)
        DECLARE_LOOKUP(14)
        DECLARE_LOOKUP(15)
        DECLARE_LOOKUP(20)
        DECLARE_LOOKUP(25)
        DECLARE_LOOKUP(30)
        DECLARE_LOOKUP(40)
        else { fprintf(stderr, "unsupported number of dimensions\n"); ASSERT2(false, "invalid dims"); }
    }
    
    void save(string filename) {
        if (0) { }
        DECLARE_SAVE(3)
        DECLARE_SAVE(4)
        DECLARE_SAVE(5)
        DECLARE_SAVE(6)
        DECLARE_SAVE(7)
        DECLARE_SAVE(8)
        DECLARE_SAVE(9)
        DECLARE_SAVE(10)
        DECLARE_SAVE(11)
        DECLARE_SAVE(12)
        DECLARE_SAVE(13)
        DECLARE_SAVE(14)
        DECLARE_SAVE(15)
        DECLARE_SAVE(20)
        DECLARE_SAVE(25)
        DECLARE_SAVE(30)
        DECLARE_SAVE(40)
        else { fprintf(stderr, "unsupported number of dimensions\n"); ASSERT2(false, "invalid dims"); }
    }
};

#endif

