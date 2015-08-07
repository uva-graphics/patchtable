
#ifndef _patch_pca_h
#define _patch_pca_h

#include "../solver/array.h"
#include <opencv2/core/core.hpp>
#include <memory>

using std::shared_ptr;

class PatchTableParams;

#define PATCH_PCA_DEFAULT_N 1000

template<class real>
void extract_patches(PatchTableParams *p, const Array<real> &a, Array<real> &patches) {
    int n = p->pca_samples;
    //printf("extract patches resize\n");
    patches.resize(n, p->patch_w*p->patch_w*a.channels());
    int aew = a.width()-p->patch_w+1;
    int aeh = a.height()-p->patch_w+1;
    //printf("extract patches begin loop %d\n", n);
    for (int i = 0; i < n; i++) {
        //printf("loop iter %d, aew, aeh=%d, %d, patches: %dx%d\n", i, aew, aeh, patches.height(), patches.width());
        int ax = rand()%aew;
        int ay = rand()%aeh;
        real *patch = &patches(i, 0);
        int count = 0;
        //printf("loop iter %d, reading patch\n", i);
        for (int dy = 0; dy < p->patch_w; dy++) {
            const real *arow = &a.get_nearest(ay+dy, ax, 0);
            for (int dx = 0; dx < p->patch_w; dx++) {
                const real *aptr = &arow[dx];
                for (int c = 0; c < a.channels(); c++) {
                    patch[count++] = aptr[c];
                }
            }
        }
        //printf("loop iter %d, done reading patch\n", i);
    }
}

template<class real>
shared_ptr<cv::PCA> get_patch_pca(PatchTableParams *p, const Array<real> &a) {
    //printf("begin get_patch_pca\n");
    Array<real> patches;
    extract_patches(p, a, patches);
    //printf("done extract_patches\n");
    
    return make_shared<cv::PCA>(patches.to_cv(), p->pca_subtract_mean ? cv::Mat(): cv::Mat::zeros(1, p->patch_w*p->patch_w*a.channels(), CV_32FC1), CV_PCA_DATA_AS_ROW, p->ndims);
}

template<class real>
void apply_patch_pca(PatchTableParams *p, const Array<real> &a, shared_ptr<cv::PCA> &pca, Array<real> &out) {
    /* Pad the same as gck.h */
    int mean_channels = 0;
    
    //printf("apply_patch_pca\n");
    out.resize(a.height()+p->patch_w-1, a.width()+p->patch_w-1, p->ndims);
    
    if (p->pca_mean) {
        Array<real> temp;
        gck<real, real>(a, temp, 1, p->nchroma > 0 ? 1: 0, p->patch_w);
        ASSERT2(temp.width() == out.width() && temp.height() == out.height(), "expected temp and out to have same size");
        mean_channels = temp.channels();
        for (int y = 0; y < temp.height(); y++) {
            for (int x = 0; x < temp.width(); x++) {
                for (int k = 0; k < temp.channels(); k++) {
                    out(y, x, k) = temp(y, x, k);
                }
            }
        }
    }

    //printf("resize arrays\n");
    cv::Mat pca_input(1, p->patch_w*p->patch_w*3, CV_32FC1);
    cv::Mat pca_result(1, p->ndims, CV_32FC1);
    float *pca_input_f = pca_input.ptr<float>(0);
    float *pca_result_f = pca_result.ptr<float>(0);
    
    //printf("write loop\n");
    for (int y = gck_ymin(out); y < gck_ymax(out); y++) {
        for (int x = gck_xmin(out); x < gck_xmax(out); x++) {
            int ay = y - gck_ymin(out);
            int ax = x - gck_xmin(out);
            //printf("write %d, %d, a size: %dx%dx%d, out size: %dx%dx%d\n", x, y, a.height(), a.width(), a.channels(), out.height(), out.width(), out.channels());
        
            /* TODO: Could optimize by performing the matrix multiply ourselves and templating on the patch width. */
            
            int count = 0;
            for (int dy = 0; dy < p->patch_w; dy++) {
                const real *arow = &a.get_nearest(ay+dy, ax, 0);
                for (int dx = 0; dx < p->patch_w; dx++) {
                    const real *aptr = &arow[dx];
                    for (int c = 0; c < a.channels(); c++) {
                        pca_input_f[count++] = aptr[c];
                    }
                }
            }
            //printf("project %d, %d\n", x, y);

            pca->project(pca_input, pca_result);
            //printf("write back result %d, %d\n", x, y);
            
            real *out_ptr = &out(y, x, 0);
            for (int i = mean_channels; i < p->ndims; i++) {
                out_ptr[i] = pca_result_f[i-mean_channels];
            }
        }
    }
    //printf("done pca\n");
}

#endif

