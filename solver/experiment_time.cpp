
#include "timer.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef __GNUC__
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

#undef ASSERT
#define ASSERT(x) if (!(x)) { fprintf(stderr, "assert failed at %s:%d\n", __FILE__, __LINE__); exit(1); }
#undef MAX
#undef MIN
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

#define TEST_W 10000
#define TEST_H 10000

bool use_random = false;

void clear_array(double *v, int n, double c=0.0) {
    for (int i = 0; i < n; i++) {
        v[i] = c;
    }
}

class Image {
public:
    double *data;
    int w, h, channels;
    
private:
    void alloc(int w_, int h_, int channels_) {
        w = w_;
        h = h_;
        channels = channels_;
        data = new double[w*h*channels];
    }

    void copy_from(const Image &other) {
        alloc(other.w, other.h, other.channels);
        memcpy((void *) data, (void *) other.data, w*h*channels*sizeof(double));
    }

public:
    Image(int w_, int h_, int channels_=3) {
        alloc(w_, h_, channels_);
    }
    
    ~Image() {
        delete[] data;
    }
    
    Image (const Image &other) {
        copy_from(other);
    }

    Image &operator = (const Image &other) {
        delete[] data;
        copy_from(other);
        return *this;
    }
    
    void clear(double c=0.0) {
        int n = w*h*channels;
        clear_array(data, n, c);
    }
    
    void random(int seed=1) {
        srand(seed);
        int n = w*h*channels;
        for (int i = 0; i < n; i++) {
            data[i] = rand() / (RAND_MAX+1.0);
        }
    }
    
    bool equals(const Image &other, double eps=0.0, int *xdiffer = NULL, int *ydiffer = NULL, int *cdiffer = NULL) {
        if (other.w != w || other.h != h || other.channels != channels) {
            return false;
        }
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < channels; c++) {
                    double a = (*this)(x,y,c);
                    double b = other(x,y,c);
                    if (fabs(a-b) >= (MAX(a, b)+1)*eps) {
                        if (xdiffer) { *xdiffer = x; }
                        if (ydiffer) { *ydiffer = y; }
                        if (cdiffer) { *cdiffer = c; }
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
#define DECLARE_INDEXING(constness) \
    ALWAYS_INLINE double &operator()(int y) constness { \
        return data[y*(w*channels)]; \
    } \
    \
    ALWAYS_INLINE double &operator()(int x, int y) constness { \
        return data[(y*w+x)*channels]; \
    } \
    \
    ALWAYS_INLINE double &operator()(int x, int y, int c) constness { \
        return data[((y*w+x)*channels)+c]; \
    }
    
    DECLARE_INDEXING();
    DECLARE_INDEXING(const);
};

#define BLUR_3_BODY() \
double *out_c = &out(x, y); \
for (int c = 0; c < out.channels; c++) { \
    out_c[c] = (in_c1[c] + 2 * in_c2[c] + in_c3[c]) / 4; \
}

#define BLUR_3_PRELUDE(step) \
double *in_c2 = &in(x, y); \
double *in_c1 = in_c2 - step; \
double *in_c3 = in_c2 + step;

#define BLUR_5_BODY() \
double *out_c = &out(x, y); \
for (int c = 0; c < out.channels; c++) { \
    out_c[c] = (in_c1[c] + in_c2[c] + in_c3[c] + in_c4[c] + in_c5[c]) / 5; \
}

#define BLUR_5_PRELUDE(step) \
double *in_c3 = &in(x, y); \
double *in_c2 = in_c3 - step; \
double *in_c1 = in_c3 - 2*step; \
double *in_c4 = in_c3 + step; \
double *in_c5 = in_c3 + 2*step; \

#define BLUR_9_BODY() \
double *out_c = &out(x, y); \
for (int c = 0; c < out.channels; c++) { \
    out_c[c] = (in_c1[c] + in_c2[c] + in_c3[c] + in_c4[c] + in_c5[c] + in_c6[c] + in_c7[c] + in_c8[c] + in_c9[c]) / 9; \
}

#define BLUR_9_PRELUDE(step) \
double *in_c5 = &in(x, y); \
double *in_c4 = in_c5 - step; \
double *in_c3 = in_c5 - 2*step; \
double *in_c2 = in_c5 - 3*step; \
double *in_c1 = in_c5 - 4*step; \
double *in_c6 = in_c5 + step; \
double *in_c7 = in_c5 + 2*step; \
double *in_c8 = in_c5 + 3*step; \
double *in_c9 = in_c5 + 4*step;

Image I_temp(1, 1);

void zero_output(Image &in, Image &out) {
    clear_array(out.data, out.w*out.h*out.channels);
}

void blur_5x5(Image &in, Image &out) {
    for (int y = 5; y < out.h-5; y++) {
        for (int x = 4; x < out.w-4; x++) {
            double *out_c = &out(x, y);
            double *row1 = &in(x,y-2);
            double *row2 = &in(x,y-1);
            double *row3 = &in(x,y);
            double *row4 = &in(x,y+1);
            double *row5 = &in(x,y+2);
            double *in_c01 = &row1[-2];
            double *in_c02 = &row1[-1];
            double *in_c03 = &row1[0];
            double *in_c04 = &row1[1];
            double *in_c05 = &row1[2];
            double *in_c06 = &row2[-2];
            double *in_c07 = &row2[-1];
            double *in_c08 = &row2[0];
            double *in_c09 = &row2[1];
            double *in_c10 = &row2[2];
            double *in_c11 = &row3[-2];
            double *in_c12 = &row3[-1];
            double *in_c13 = &row3[0];
            double *in_c14 = &row3[1];
            double *in_c15 = &row3[2];
            double *in_c16 = &row4[-2];
            double *in_c17 = &row4[-1];
            double *in_c18 = &row4[0];
            double *in_c19 = &row4[1];
            double *in_c20 = &row4[2];
            double *in_c21 = &row5[-2];
            double *in_c22 = &row5[-1];
            double *in_c23 = &row5[0];
            double *in_c24 = &row5[1];
            double *in_c25 = &row5[2];
            for (int c = 0; c < out.channels; c++) {
                out_c[c] = (1.0/25.0) * (in_c01[c] +
                                         in_c02[c] +
                                         in_c03[c] +
                                         in_c04[c] +
                                         in_c05[c] +
                                         in_c06[c] +
                                         in_c07[c] +
                                         in_c08[c] +
                                         in_c09[c] +
                                         in_c10[c] +
                                         in_c11[c] +
                                         in_c12[c] +
                                         in_c13[c] +
                                         in_c14[c] +
                                         in_c15[c] +
                                         in_c16[c] +
                                         in_c17[c] +
                                         in_c18[c] +
                                         in_c19[c] +
                                         in_c20[c] +
                                         in_c21[c] +
                                         in_c22[c] +
                                         in_c23[c] +
                                         in_c24[c] +
                                         in_c25[c]);
            }
        }
    }
}

void blur_25_horiz(Image &in, Image &out) {
    for (int y = 0; y < out.h; y++) {
        for (int x = 12; x < out.w-12; x++) {
            double *out_c = &out(x, y);
            for (int c = 0; c < out.channels; c++) {
                double accum = 0.0;
                for (int dx = -12; dx <= 12; dx++) {
                    accum += in(x+dx,y,c);
                }
                out_c[c] = (1.0/25.0) * accum;
            }
        }
    }
}

void blur_25_vert(Image &in, Image &out) {
    for (int y = 12; y < out.h-12; y++) {
        for (int x = 0; x < out.w; x++) {
            double *out_c = &out(x, y);
            for (int c = 0; c < out.channels; c++) {
                double *row0 = &in(x,y,c);
                int stride = out.w*out.channels;
                double accum = 0.0;
                for (int dy = -12; dy <= 12; dy++) {
                    accum += row0[dy*stride];
                }
                out_c[c] = (1.0/25.0) * accum;
            }
        }
    }
}

void blur_5x5_sparse(Image &in, Image &out) {
    for (int y = 5; y < out.h-5; y++) {
        for (int x = 4; x < out.w-4; x++) {
            double *out_c = &out(x, y);
            double *row1 = &in(x,y-2);
            double *row2 = &in(x,y-1);
            double *row3 = &in(x,y);
            double *row4 = &in(x,y+1);
            double *row5 = &in(x,y+2);
            double *in_c01 = &row1[-4];
            double *in_c02 = &row1[-2];
            double *in_c03 = &row1[0];
            double *in_c04 = &row1[2];
            double *in_c05 = &row1[4];
            double *in_c06 = &row2[-4];
            double *in_c07 = &row2[-2];
            double *in_c08 = &row2[0];
            double *in_c09 = &row2[2];
            double *in_c10 = &row2[4];
            double *in_c11 = &row3[-4];
            double *in_c12 = &row3[-2];
            double *in_c13 = &row3[0];
            double *in_c14 = &row3[2];
            double *in_c15 = &row3[4];
            double *in_c16 = &row4[-4];
            double *in_c17 = &row4[-2];
            double *in_c18 = &row4[0];
            double *in_c19 = &row4[2];
            double *in_c20 = &row4[4];
            double *in_c21 = &row5[-4];
            double *in_c22 = &row5[-2];
            double *in_c23 = &row5[0];
            double *in_c24 = &row5[2];
            double *in_c25 = &row5[4];
            for (int c = 0; c < out.channels; c++) {
                out_c[c] = (1.0/25.0) * (in_c01[c] +
                                         in_c02[c] +
                                         in_c03[c] +
                                         in_c04[c] +
                                         in_c05[c] +
                                         in_c06[c] +
                                         in_c07[c] +
                                         in_c08[c] +
                                         in_c09[c] +
                                         in_c10[c] +
                                         in_c11[c] +
                                         in_c12[c] +
                                         in_c13[c] +
                                         in_c14[c] +
                                         in_c15[c] +
                                         in_c16[c] +
                                         in_c17[c] +
                                         in_c18[c] +
                                         in_c19[c] +
                                         in_c20[c] +
                                         in_c21[c] +
                                         in_c22[c] +
                                         in_c23[c] +
                                         in_c24[c] +
                                         in_c25[c]);
            }
        }
    }
}

void blur_5x5_sparse_b(Image &in, Image &out) {
    for (int y = 5; y < out.h-5; y++) {
        for (int x = 4; x < out.w-4; x++) {
            double *out_c = &out(x, y);
            double *row1 = &in(x,y-4);
            double *row2 = &in(x,y-3);
            double *row3 = &in(x,y-2);
            double *row4 = &in(x,y-1);
            double *row5 = &in(x,y);
            double *row6 = &in(x,y+1);
            double *row7 = &in(x,y+2);
            double *row8 = &in(x,y+3);
            double *row9 = &in(x,y+4);
            double *row10 = &in(x,y+5);
            
            double *in_c01 = &row1[-4];
            double *in_c02 = &row2[-2];
            double *in_c03 = &row1[0];
            double *in_c04 = &row2[2];
            double *in_c05 = &row1[4];
            
            double *in_c06 = &row3[-4];
            double *in_c07 = &row4[-2];
            double *in_c08 = &row3[0];
            double *in_c09 = &row4[2];
            double *in_c10 = &row3[4];
            
            double *in_c11 = &row5[-4];
            double *in_c12 = &row6[-2];
            double *in_c13 = &row5[0];
            double *in_c14 = &row6[2];
            double *in_c15 = &row5[4];
            
            double *in_c16 = &row7[-4];
            double *in_c17 = &row8[-2];
            double *in_c18 = &row7[0];
            double *in_c19 = &row8[2];
            double *in_c20 = &row7[4];
            
            double *in_c21 = &row9[-4];
            double *in_c22 = &row10[-2];
            double *in_c23 = &row9[0];
            double *in_c24 = &row10[2];
            double *in_c25 = &row9[4];
            
            for (int c = 0; c < out.channels; c++) {
                out_c[c] = (1.0/25.0) * (in_c01[c] +
                                         in_c02[c] +
                                         in_c03[c] +
                                         in_c04[c] +
                                         in_c05[c] +
                                         in_c06[c] +
                                         in_c07[c] +
                                         in_c08[c] +
                                         in_c09[c] +
                                         in_c10[c] +
                                         in_c11[c] +
                                         in_c12[c] +
                                         in_c13[c] +
                                         in_c14[c] +
                                         in_c15[c] +
                                         in_c16[c] +
                                         in_c17[c] +
                                         in_c18[c] +
                                         in_c19[c] +
                                         in_c20[c] +
                                         in_c21[c] +
                                         in_c22[c] +
                                         in_c23[c] +
                                         in_c24[c] +
                                         in_c25[c]);
            }
        }
    }
}

void blur_3_horiz(Image &in, Image &out) {
    for (int y = 0; y < out.h; y++) {
        for (int x = 1; x < out.w-1; x++) {
            BLUR_3_PRELUDE(out.channels);
            BLUR_3_BODY();
        }
    }
}

void blur_3_vert(Image &in, Image &out) {
    int stride = out.w*out.channels;
    for (int y = 1; y < out.h-1; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_3_PRELUDE(stride);
            BLUR_3_BODY();
        }
    }
}

void blur_5_horiz(Image &in, Image &out) {
    for (int y = 0; y < out.h; y++) {
        for (int x = 2; x < out.w-2; x++) {
            BLUR_5_PRELUDE(out.channels);
            BLUR_5_BODY();
        }
    }
}

void blur_5_vert(Image &in, Image &out) {
    int stride = out.w*out.channels;
    for (int y = 2; y < out.h-2; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_5_PRELUDE(stride);
            BLUR_5_BODY();
        }
    }
}

void blur_9_horiz(Image &in, Image &out) {
    for (int y = 0; y < out.h; y++) {
        for (int x = 4; x < out.w-4; x++) {
            BLUR_9_PRELUDE(out.channels);
            BLUR_9_BODY();
        }
    }
}

void blur_9_vert(Image &in, Image &out) {
    int stride = out.w*out.channels;
    for (int y = 4; y < out.h-4; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_9_PRELUDE(stride);
            BLUR_9_BODY();
        }
    }
}

void blur_3_separable(Image &in, Image &out) {
    blur_3_horiz(in, I_temp);
    blur_3_vert(I_temp, out);
}

void blur_5_separable(Image &in, Image &out) {
    blur_5_horiz(in, I_temp);
    blur_5_vert(I_temp, out);
}

void blur_9_separable(Image &in, Image &out) {
    blur_9_horiz(in, I_temp);
    blur_9_vert(I_temp, out);
}

void blur_3_sparse_horiz(Image &in, Image &out) {
    int step = 10;
    for (int y = 0; y < out.h; y++) {
        for (int x = step; x < out.w-step; x++) {
            BLUR_3_PRELUDE(out.channels*step);
            BLUR_3_BODY();
        }
    }
}

void blur_3_sparse_vert(Image &in, Image &out) {
    int step = 10;
    int stride = out.w*out.channels*step;
    for (int y = step; y < out.h-step; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_3_PRELUDE(stride);
            BLUR_3_BODY();
        }
    }
}

void blur_5_sparse_horiz(Image &in, Image &out) {
    int step = 5;
    for (int y = 0; y < out.h; y++) {
        for (int x = step*2; x < out.w-step*2; x++) {
            BLUR_5_PRELUDE(out.channels*step);
            BLUR_5_BODY();
        }
    }
}

void blur_5_sparse_vert(Image &in, Image &out) {
    int step = 5;
    int stride = out.w*out.channels*step;
    for (int y = step*2; y < out.h-step*2; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_5_PRELUDE(stride);
            BLUR_5_BODY();
        }
    }
}

void blur_9_sparse_horiz(Image &in, Image &out) {
    int step = 5;
    for (int y = 0; y < out.h; y++) {
        for (int x = step*4; x < out.w-step*4; x++) {
            BLUR_9_PRELUDE(out.channels*step);
            BLUR_9_BODY();
        }
    }
}

void blur_9_sparse_vert(Image &in, Image &out) {
    int step = 5;
    int stride = out.w*out.channels*step;
    for (int y = step*4; y < out.h-step*4; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_9_PRELUDE(stride);
            BLUR_9_BODY();
        }
    }
}

#define BLUR_IIR_BODY(PREV_OUT) \
double *in_c = &in(x, y); \
double *out_c = &out(x, y); \
for (int channel = 0; channel < out.channels; channel++) { \
    accum[channel] = accum[channel]*0.5 + in_c[channel]; \
    out_c[channel] = PREV_OUT + accum[channel]; \
}

void blur_2_iir_horiz(Image &in, Image &out) {
    double *accum = new double[out.channels];
    
    for (int y = 0; y < out.h; y++) {
        clear_array(accum, out.channels);
        for (int x = 0; x < out.w; x++) {
            BLUR_IIR_BODY(0.0);
        }
        clear_array(accum, out.channels);
        for (int x = out.w-1; x >= 0; x--) {
            BLUR_IIR_BODY(out_c[channel]);
        }
    }
    delete[] accum;
}

void blur_2_iir_vert(Image &in, Image &out) {
    double *accum = new double[out.channels];
    
    for (int x = 0; x < out.w; x++) {
        clear_array(accum, out.channels);
        for (int y = 0; y < out.h; y++) {
            BLUR_IIR_BODY(0.0);
        }
        clear_array(accum, out.channels);
        for (int y = out.h-1; y >= 0; y--) {
            BLUR_IIR_BODY(out_c[channel]);
        }
    }
    delete[] accum;
}

#define BLUR_IIR_BODY2(PREV_OUT) \
double *in_c = &in(x, y); \
double *out_c = &out(x, y); \
double *accum = &accum_row[x*out.channels]; \
for (int channel = 0; channel < out.channels; channel++) { \
    accum[channel] = accum[channel]*0.5 + in_c[channel]; \
    out_c[channel] = PREV_OUT + accum[channel]; \
}

void blur_2_iir_vert2(Image &in, Image &out) {
    int nrow = out.channels*out.w;
    double *accum_row = new double[nrow];
    
    clear_array(accum_row, nrow);
    for (int y = 0; y < out.h; y++) {
        for (int x = 0; x < out.w; x++) {
            BLUR_IIR_BODY2(0.0);
        }
    }
    clear_array(accum_row, nrow);
    for (int y = out.h-1; y >= 0; y--) {
        for (int x = 0; x < out.w; x++) {
            BLUR_IIR_BODY2(out_c[channel]);
        }
    }
    delete[] accum_row;
}

typedef void (*FuncType)(Image &, Image &);

void init_test_images(Image &I, Image &Iout) {
    if (use_random) {
        I.random();
    } else {
        I.clear(0.5);
    }
    Iout.clear(-1.0);
    I(0,1,0) = 1.0;
    I(2,1,0) = 0.6;
    I(1,0,0) = 1.0;
    I(1,2,0) = 0.7;
}

double test_func(FuncType f) {
    int w = TEST_W;
    int h = TEST_H;
    Image I(w, h);
    Image Iout(w, h);
    I_temp = I;
    init_test_images(I, Iout);
    
    double T0 = wall_time();
    f(I, Iout);
    double T = wall_time()-T0;
    return T;
}

void test_iir() {
    int w = TEST_W;
    int h = TEST_H;
    Image I(w, h);
    Image Iout(w, h);
    Image Iout2(w, h);
    init_test_images(I, Iout);
    init_test_images(I, Iout2);
    
    blur_2_iir_vert(I, Iout);
    blur_2_iir_vert2(I, Iout2);
    
    int xdiffer = 0, ydiffer = 0, cdiffer = 0;
    bool is_equal = Iout.equals(Iout2, 1e-8, &xdiffer, &ydiffer, &cdiffer);
    if (!is_equal) {
        printf("Images differ\n");
        printf("  at %d,%d,%d\n", xdiffer, ydiffer, cdiffer);
        printf("  %.16lf %.16lf\n", Iout(xdiffer, ydiffer, cdiffer), Iout2(xdiffer, ydiffer, cdiffer));
    }
}

const char *func_names[] = {
    "blur_5x5",
    "blur_5x5_sparse",
    "blur_5x5_sparse_b",
    "blur_25_horiz",
    "blur_25_vert",
    "zero_output",
    "blur_3_horiz",
    "blur_3_vert",
    "blur_5_horiz",
    "blur_5_vert",
    "blur_9_horiz",
    "blur_9_vert",
    "blur_3_sparse_horiz",
    "blur_3_sparse_vert",
    "blur_5_sparse_horiz",
    "blur_5_sparse_vert",
    "blur_9_sparse_horiz",
    "blur_9_sparse_vert",
    "blur_2_iir_horiz",
    "blur_2_iir_vert2",
    "blur_3_separable",
    "blur_5_separable",
    "blur_9_separable",
};

FuncType func_list[] = {
    blur_5x5,
    blur_5x5_sparse,
    blur_5x5_sparse_b,
    blur_25_horiz,
    blur_25_vert,
    zero_output,
    blur_3_horiz,
    blur_3_vert,
    blur_5_horiz,
    blur_5_vert,
    blur_9_horiz,
    blur_9_vert,
    blur_3_sparse_horiz,
    blur_3_sparse_vert,
    blur_5_sparse_horiz,
    blur_5_sparse_vert,
    blur_9_sparse_horiz,
    blur_9_sparse_vert,
    blur_2_iir_horiz,
    blur_2_iir_vert2,
    blur_3_separable,
    blur_5_separable,
    blur_9_separable,
};

int func_taps[] = {         /* Total taps used in all passes */
    25,
    25,
    25,
    25,
    25,
    0,
    3,
    3,
    5,
    5,
    9,
    9,
    3,
    3,
    5,
    5,
    9,
    9,
    4,
    4,
    6,
    10,
    18,
};

int func_passes[] = {
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
};

int func_surfacex[] = {
    10,
    30,
    40,
    2,
    50,
    0,
    2,
    6,
    2,
    10,
    2,
    18,
    6,
    6,
    10,
    10,
    18,
    18,
    4,
    6,
    8,
    12,
    20,
};

int func_surfacey[] = {
    10,
    10,
    40,
    50,
    2,
    0,
    6,
    2,
    10,
    2,
    18,
    2,
    6,
    6,
    10,
    10,
    18,
    18,
    6,
    4,
    8,
    12,
    20,
};

int nfuncs = 23;

int main() {
//    test_iir();
    
    for (int i = 0; i < nfuncs; i++) {
        printf("%-20s %d %d %d %d %f\n", func_names[i], func_taps[i], func_passes[i], func_surfacex[i], func_surfacey[i], test_func(func_list[i]));
    }
    
    return 0;
}
