
/* Templated arbitrary-dimensional rectangular array with efficient lookup. */

#ifndef _array_h
#define _array_h

#define cimg_display 0

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string.h>
#include <assert.h>
//#include "CImg.h"
#include "util.h"
#include "timer.h"
#include "pfm.h"
#include <png++/png.hpp>

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/map.hpp> 
#include <boost/mpl/at.hpp>

#include <unordered_set>

#define USE_HALIDE 0
#define USE_ADEPT  0
#ifndef ARRAY_USE_OPENCV
#define ARRAY_USE_OPENCV 0
#endif

#if ARRAY_USE_OPENCV
#include <opencv2/opencv.hpp>
#endif
   
#define DEBUG_TIME         0                    /* Print running time information about various functions */
#define DEBUG_SAVE_IMAGES  0                    /* Save some intermediate images */

#if USE_HALIDE
#include <Halide.h>
using Halide::Image;
#include "../../Halide/apps/support/image_io.h"
#endif

#if USE_ADEPT
#include <adept.h>
using adept::adouble;
#endif

using std::string;
using std::vector;
using std::unordered_set;
using std::set;
using std::pair;
using std::map;
//using cimg_library::CImg;

//extern bool array_debug;

#ifdef __GNUC__
#define INLINE inline __attribute__((always_inline))
#else
#define INLINE __forceinline
#endif

#undef MIN
#undef MAX
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))
#define in_bounds(a, b) ((unsigned) (a) < (unsigned) (b))
#define CLAMP_INT(x) (MAX(0,MIN(x,255)))

#define DEBUG BUILD_DEBUG
#define ASSERT(x, msg) if (DEBUG && !(x)) { fprintf(stderr, "Assertion failed in %s(%d): %s\n", __FILE__, __LINE__, msg); assert(false); exit(1); }
#define rand_f() (rand()/(RAND_MAX+1.0))
#define rand_uniform() (float(rand())/float(RAND_MAX+1.0))

#define R_COORD 0
#define G_COORD 1
#define B_COORD 2 

#define X_COORD 0
#define Y_COORD 1
#define Z_COORD 2 

#define L_COORD 0
#define A_COORD 1
#define B_COORD 2 

#if ARRAY_USE_OPENCV

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,CV_64FC1> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,CV_32FC1> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<int,CV_8UC1> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<int,CV_16UC1> >
  > cvTypeChannels1;

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,CV_64FC2> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,CV_32FC2> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<int,CV_8UC2> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<int,CV_16UC2> >
  > cvTypeChannels2;

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,CV_64FC3> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,CV_32FC3> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<int,CV_8UC3> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<int,CV_16UC3> >
  > cvTypeChannels3;

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,CV_64FC4> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,CV_32FC4> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<int,CV_8UC4> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<int,CV_16UC4> >
  > cvTypeChannels4;

#define UNKNOWN_CV_TYPE -10000

#define CV_TYPE(I, real) ((I).channels() == 1 ? boost::mpl::at<cvTypeChannels1, real>::type::value: \
                         ((I).channels() == 2 ? boost::mpl::at<cvTypeChannels2, real>::type::value: \
                         ((I).channels() == 3 ? boost::mpl::at<cvTypeChannels3, real>::type::value: \
                         ((I).channels() == 4 ? boost::mpl::at<cvTypeChannels4, real>::type::value: UNKNOWN_CV_TYPE) \
                         )))

#endif

template<class real>
string vector_to_str_real(const vector<real> &v) {
    string ans("[");
    char buf[256];
    for (int i = 0; i < (int) v.size(); i++) {
        sprintf(buf, "%e", double(v[i]));
        ans += buf;
        if (i < (int) v.size() - 1) {
            ans += ",";
        }
    }
    ans += "]";
    return ans;
}

template<class int_t>
string vector_to_str_int(const vector<int_t> &v) {
    string ans("[");
    char buf[256];
    for (int i = 0; i < (int) v.size(); i++) {
        sprintf(buf, "%d", v[i]);
        ans += buf;
        if (i < (int) v.size() - 1) {
            ans += ",";
        }
    }
    ans += "]";
    return ans;
}

template<class U, class V>
INLINE V cast_types(U x) {
    return x;
}

#if USE_ADEPT
template<>
INLINE adouble cast_types(adouble x) {
    return x;
}

template<>
INLINE double cast_types(adouble x) {
    return x.value();
}

template<>
INLINE adouble cast_types(double x) {
    return x;
}

INLINE double to_double(adouble x) {
    return x.value();
}
#endif

template<>
INLINE double cast_types(double x) {
    return x;
}

INLINE double to_double(double x) {
    return x;
}

//INLINE double abs(double x) {
//    return fabs(x);
//}

INLINE string num_to_str(double x) {
    char buf[256];
    sprintf(buf, "%11.6g", x);
    return buf;
}

INLINE string num_to_str(int x) {
    char buf[256];
    sprintf(buf, "%d", x);
    return string(buf);
}

#if USE_ADEPT
INLINE string num_to_str(adouble x) {
    return num_to_str(x.value());
}
#endif

template<class real, int N=1024>
class PowTable { public:
    real table[N];
    real upper;
    real scale;
    PowTable(real p, real upper_=1) :upper(upper_) {
        scale = N/upper;
        for (int i = 0; i < N; i++) {
            double v = (i+0.5)/N*upper;
            table[i] = pow(v, p);
        }
    }
    real operator() (real v) {
        int i = int(v*scale);
        if (i < 0) { i = 0; }
        else if (i >= N) { i = N-1; }
        return table[i];
    }
};

#define DIV_ROUND_UP(x, y) (((x)+(y)-1)/(y))
#if VISUAL_STUDIO_WORKAROUND
#define TEMPLATE_OTHER_ARRAY template<class T>
#define OTHER_ARRAY Array<T>
#else
#define TEMPLATE_OTHER_ARRAY template<class T, int T_vectorize=-1, int T_downsample_dim0=1, int T_downsample_dim1=1, int T_downsample_dim2=1>
#define OTHER_ARRAY Array<T, T_vectorize, T_downsample_dim0, T_downsample_dim1, T_downsample_dim2>
#endif

template<class real, int vectorize=-1, int downsample_dim0=1, int downsample_dim1=1, int downsample_dim2=1, int downsample_dim3=1>
class Array {
    public:
    
    real *data;
    vector<int> sizes;                  /* Sizes at full resolution */
    vector<int> sizes_downsampled;      /* Downsampled sizes */
    vector<int> stride;
    int nelems;
    bool dealloc;
    
    /* Changes array sizes. Any existing data is lost -- new data will be uninitialized. */
    void resize(const vector<int> &sizes_, void *data_ptr=NULL) {
        /* TODO: If downsampled allocate less memory, and decrease nelems accordingly. */
        if (sizes != sizes_) {
            if (vectorize >= 0) {
                ASSERT(vectorize < sizes_.size(), "expected vectorize to be < # number of dims");
            }
//            printf("Array::resize %s %s\n", vector_to_str_int(sizes).c_str(), vector_to_str_int(sizes_).c_str());
//            if (array_debug) { ASSERT2(false, "array should not be resized here\n"); }
            if (dealloc) {
                delete[] data;
            }
            sizes = sizes_;
            
            sizes_downsampled = sizes_;
            if (sizes_downsampled.size() > 0 && downsample_dim0 != 1) { sizes_downsampled[0] = DIV_ROUND_UP(sizes_downsampled[0], downsample_dim0); }
            if (sizes_downsampled.size() > 1 && downsample_dim1 != 1) { sizes_downsampled[1] = DIV_ROUND_UP(sizes_downsampled[1], downsample_dim1); }
            if (sizes_downsampled.size() > 2 && downsample_dim2 != 1) { sizes_downsampled[2] = DIV_ROUND_UP(sizes_downsampled[2], downsample_dim2); }
            
            stride.resize(sizes.size());
            nelems = 1;
            for (int i = stride.size()-1; i >= 0; i--) {
                stride[i] = nelems;
                nelems *= (i != vectorize) ? sizes[i]: 1;
            }
            if (!data_ptr) {
                data = new real[nelems];
                dealloc = true;
            } else {
                data = (real *) data_ptr;
                dealloc = false;
            }
        }
    }

    void resize(int d0, int d1, int d2, int d3) {
        vector<int> L;
        L.push_back(d0);
        L.push_back(d1);
        L.push_back(d2);
        L.push_back(d3);
        resize(L);
    }
    
    void resize(int h, int w, int depth) {
        vector<int> L;
        L.push_back(h);
        L.push_back(w);
        L.push_back(depth);
        resize(L);
    }

    void resize(int h, int w) {
        vector<int> L;
        L.push_back(h);
        L.push_back(w);
        resize(L);
    }

    void resize(int h) {
        vector<int> L;
        L.push_back(h);
        resize(L);
    }

    vector<real> as_vector() const {
        ASSERT(dimensions() == 1, "as_vector expected 1 dimensional array");
        int n = size();
        vector<real> ans(n);
        for (int i = 0; i < n; i++) {
            ans[i] = (*this)(i);
        }
        return ans;
    }
    
    void resize(const Array<int> &sizes_) {
        resize(sizes_.as_vector());
    }
    
    Array() :data(NULL), dealloc(true) {
		vector<int> sizeL;
		sizeL.push_back(1);
        if (vectorize >= 0) {
            while (vectorize >= sizeL.size()) {
                sizeL.push_back(1);
            }
        }
        resize(sizeL);
    }
    
//    Array(const CImg<real> &image) :data(NULL), dealloc(true) { // NOTE: This only grabs the red channel
//    	resize({ image.height(), image.width() });
//		for (int y = 0; y < image.height() ; y++){
//			for (int x = 0; x < image.width() ; x++){
//				(*this)(y,x) = image(x,y,0,0);
//			}
//		}
//    }
    
#if USE_HALIDE
    Array(const Halide::Image<real> &image) :data(NULL), dealloc(true) {
    	resize({ image.height(), image.width(), image.channels() });
	    for (int z = 0; z < image.channels() ; z++){
		    for (int y = 0; y < image.height() ; y++){
			    for (int x = 0; x < image.width() ; x++){
				    (*this)(y,x,z) = image(x,y,z);
			    }
		    }
	    }
    }
#endif

    Array(int h_, int w_, int depth_) :data(NULL), dealloc(true) {
        resize(h_, w_, depth_);
    }

    Array(int h_, int w_) :data(NULL), dealloc(true) {
        resize(h_, w_);
    }

    Array(int h_) :data(NULL), dealloc(true) {
        resize(h_);
    }

    Array(const vector<int> &sizes_, void *data_ptr = NULL) :data(NULL), dealloc(true) {
        resize(sizes_, data_ptr);
    }

    ~Array() {
        if (dealloc) {
            delete[] data;
        }
    }

    /**
     * Writes the image to the file specified by `filename`
     * in a human-readible, ASCII format as follows:
     *
     *  - The first line is a header consisting of 3 integers specifying the
     *    width, height, and channels of the image, respectively.
     *
     *  - channels() * height() lines, each with width() space-delimited
     *    floats, formatted according to std::iostream.
     *
     * See loadAscii() for reading these files.
     */
    void saveAscii(string filename) const {
        std::ofstream f;

        f.open(filename.c_str());

        int dim = sizes.size();

        int width, height, channels;

        height   = (dim >= 1) ? this->height()   : 1;
        width    = (dim >= 2) ? this->width()    : 1;
        channels = (dim >= 3) ? this->channels() : 1;

        f << width << " ";
        f << height << " ";
        f << channels << " ";

        f << std::endl;

        for (int c = 0; c < channels; c++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    float val;

                    if (dim == 1) {
                        val = (*this)(y);
                    } else if (dim == 2) {
                        val = (*this)(y, x);
                    } else if (dim == 3) {
                        val = (*this)(y, x, c);
                    }

                    f << val << " ";
                }

                f << std::endl;
            }

            f << std::endl;
        }

        f.close();
    }

    /**
     * Reads an image from an ASCII-format file written via saveAscii().
     */
    void loadAscii(string filename) {
        FILE *f = fopen(filename.c_str(), "rt");

        int width = 0, height = 0, depth = 0;
        if (fscanf(f, "%d %d %d", &width, &height, &depth) != 3) {
            fprintf(stderr, "could not read width height depth from %s\n", filename.c_str()); exit(1);
        }
        
        resize(height, width, depth);

        int dims = sizes.size();

        for (int c = 0; c < depth; c++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    float val = 0.0;

                    if (fscanf(f, "%f", &val) != 1) {
                        fprintf(stderr, "could not read %d x %d x %d floats from %s\n", width, height, depth, filename.c_str()); exit(1);
                    }

                    if (dims == 3) {
                        get_nearest(y, x, c) = val;
                    } else {
                        get_nearest(y, x) = val;
                    }
                }
            }
        }

        fclose(f);
    }
    
    TEMPLATE_OTHER_ARRAY
    void assign(const OTHER_ARRAY &other) {
        resize(other.sizes);
#if !VISUAL_STUDIO_WORKAROUND
        if (downsample_dim0 == T_downsample_dim0 && downsample_dim1 == T_downsample_dim1 && downsample_dim2 == T_downsample_dim2) {
#endif
            for (int i = 0; i < nelems; i++) {
                data[i] = cast_types<T, real>(other.data[i]);
            }
#if !VISUAL_STUDIO_WORKAROUND
        } else {
            int dims = dimensions();
            if (dims == 1) {
                for (int v0 = 0; v0 < sizes[0]; v0 += downsample_dim0) {
                    get_nearest(v0) = other(v0);
                }
            } else if (dims == 2) {
                for (int v0 = 0; v0 < sizes[0]; v0 += downsample_dim0) {
                    for (int v1 = 0; v1 < sizes[1]; v1 += downsample_dim1) {
                        get_nearest(v0, v1) = other(v0, v1);
                    }
                }
            } else if (dims == 3) {
                for (int v0 = 0; v0 < sizes[0]; v0 += downsample_dim0) {
                    for (int v1 = 0; v1 < sizes[1]; v1 += downsample_dim1) {
                        for (int v2 = 0; v2 < sizes[2]; v2 += downsample_dim2) {
                            get_nearest(v0, v1, v2) = other(v0, v1, v2);
                        }
                    }
                }
            } else {
                ASSERT2(false, "dims > 3 not implemented");
            }
        }
#endif
    }

    Array &operator=(const Array &other) {
        assign(other);
        return *this;
    }
    
    template<class T_other>
    Array &operator =(const T_other &other) {
        assign(other);
        return *this;
    }

    Array &operator =(const vector<int> &sizes_) {
        resize(sizes_);
        return *this;
    }

    TEMPLATE_OTHER_ARRAY
    Array(const OTHER_ARRAY &other) :data(NULL), dealloc(true) {
        assign(other);
    }

    template<class Array_Other>
    Array(const Array_Other &other) :data(NULL), dealloc(true) {
        assign(other);
    }
    
    Array(const Array<real> &other) :data(NULL), dealloc(true) {
        assign(other);
    }
    
    INLINE int size() const {
        return sizes[0];
    }
    
    INLINE int height() const {
        ASSERT(sizes.size() >= 2, "need at least 2 dims to get height");
        return sizes[0];
    }

    INLINE int width() const {
        ASSERT(sizes.size() >= 2, "need at least 2 dims to get height");
        return sizes[1];
    }

    INLINE int channels() const {
        //ASSERT(sizes.size() >= 3, "need at least 3 dims to get channels");
        return (sizes.size() >= 3) ? sizes[2] : 1;
    }
    
    INLINE int dimensions() const {
        return sizes.size();
    }
    
    vector<int> get_sizes() const {                      /* Safer version of sizes, returning 1D vector<int> */
        return sizes;
    }
    
    /* get_unsafe(y, x) does nearest neighbor lookup without bounds check */
    
    INLINE real &get_unsafe(int v0) const {
        v0 /= downsample_dim0;
        ASSERT(sizes.size() == 1, "1D lookup in array of other dimensionality");
        return data[v0];
    }
    
    INLINE real &get_unsafe(int v0, int v1) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        ASSERT(sizes.size() == 2, "2D lookup in array of other dimensionality");
        return data[v0*stride[0]+v1];
    }
    
    INLINE real &get_unsafe(int v0, int v1, int v2) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        v2 /= downsample_dim2;
        ASSERT(sizes.size() == 3, "3D lookup in array of other dimensionality");
        return data[v0*stride[0]+v1*stride[1]+v2];
    }

    INLINE real &get_unsafe(int v0, int v1, int v2, int v3) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        v2 /= downsample_dim2;
        v3 /= downsample_dim3;
        ASSERT(sizes.size() == 4, "4D lookup in array of other dimensionality");
        return data[v0*stride[0]+v1*stride[1]+v2*stride[2]+v3];
    }

    /* get_nearest(y, x) does nearest neighbor lookup with bounds check in debug mode */

#define PRINT_COORDS(v0, v1) (v0 == 2 && v1 == 2)
    
    INLINE real &get_nearest(int v0) const {
        v0 /= downsample_dim0;
        ASSERT(sizes.size() >= 1, "1D lookup in array of lesser dimensionality");
        ASSERT(v0 >= 0 && v0 < sizes_downsampled[0], "1D lookup out of bounds");
        return data[v0];
    }
    
    INLINE real &get_nearest(int v0, int v1) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        ASSERT(sizes.size() >= 2, "2D lookup in array of lesser dimensionality");
        ASSERT(v0 >= 0 && v0 < sizes_downsampled[0] && v1 >= 0 && v1 < sizes_downsampled[1], "2D lookup out of bounds");
        return data[v0*stride[0]+v1];
    }
    
    INLINE real &get_nearest(int v0, int v1, int v2) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        v2 /= downsample_dim2;
        ASSERT(sizes.size() >= 3, "3D lookup in array of lesser dimensionality");
        ASSERT(v0 >= 0 && v0 < sizes_downsampled[0] && v1 >= 0 && v1 < sizes_downsampled[1] && v2 >= 0 && v2 < sizes_downsampled[2], "3D lookup out of bounds");
        return data[v0*stride[0]+v1*stride[1]+v2];
    }

    INLINE real &get_nearest(int v0, int v1, int v2, int v3) const {
        v0 /= downsample_dim0;
        v1 /= downsample_dim1;
        v2 /= downsample_dim2;
        v3 /= downsample_dim3;
        ASSERT(sizes.size() >= 4, "4D lookup in array of lesser dimensionality");
        ASSERT(v0 >= 0 && v0 < sizes_downsampled[0] && v1 >= 0 && v1 < sizes_downsampled[1] && v2 >= 0 && v2 < sizes_downsampled[2] && v3 >= 0 && v3 < sizes_downsampled[3], "4D lookup out of bounds");
        return data[v0*stride[0]+v1*stride[1]+v2*stride[2]+v3];
    }

#define ARRAY_INDEX_TO_INT() \
        int lookup = 0; \
        for (int i = 0; i < (int) index.size(); i++) { \
            ASSERT(in_bounds(index[i], sizes_downsampled[i]), "nd-lookup out of bounds"); \
            lookup += index[i]*stride[i]; \
        }

    INLINE real &get_nearest(const vector<int> &index) const {
        ASSERT(sizes.size() == index.size(), "nd-lookup in array of other dimensionality");
        ARRAY_INDEX_TO_INT();
        return data[lookup];
    }
    
    INLINE int index_to_int(const vector<int> &index) const {
        ARRAY_INDEX_TO_INT();
        return lookup;
    }
    
#include "array_copy_rect.h"

    #define CLAMP_COORD(v, dim) if ((unsigned) v >= (unsigned) sizes[dim]) { if (v < 0) { v = 0; } else if (v >= sizes[dim]) { v = sizes[dim]-1; } }
    
    INLINE real get_zero(int v0) const {
        if (in_bounds(v0, sizes_downsampled[0])) {
            return get_nearest(v0);
        }
        return get_nearest(0)*0;
    }

    INLINE real get_zero(int v0, int v1) const {
        if (in_bounds(v0, sizes_downsampled[0]) && in_bounds(v1, sizes_downsampled[1])) {
            return get_nearest(v0, v1);
        }
        return get_nearest(0, 0)*0;
    }

    INLINE real get_zero(int v0, int v1, int v2) const {
        if (in_bounds(v0, sizes_downsampled[0]) && in_bounds(v1, sizes_downsampled[1]) && in_bounds(v2, sizes_downsampled[2])) {
            return get_nearest(v0, v1, v2);
        }
        return get_nearest(0, 0, 0)*0;
    }

    INLINE real get_zero(int v0, int v1, int v2, int v3) const {
        if (in_bounds(v0, sizes_downsampled[0]) && in_bounds(v1, sizes_downsampled[1]) && in_bounds(v2, sizes_downsampled[2]) && in_bounds(v3, sizes_downsampled[3])) {
            return get_nearest(v0, v1, v2, v3);
        }
        return get_nearest(0, 0, 0, 0)*0;
    }
    
    INLINE real &get_clamp(int v0) const {
        CLAMP_COORD(v0, 0);
        return get_nearest(v0);
    }

    INLINE real &get_clamp(int v0, int v1) const {
        CLAMP_COORD(v0, 0);
        CLAMP_COORD(v1, 1);
        return get_nearest(v0, v1);
    }

    INLINE real &get_clamp(int v0, int v1, int v2) const {
        CLAMP_COORD(v0, 0);
        CLAMP_COORD(v1, 1);
        CLAMP_COORD(v2, 2);
        return get_nearest(v0, v1, v2);
    }

    INLINE real &get_clamp(int v0, int v1, int v2, int v3) const {
        CLAMP_COORD(v0, 0);
        CLAMP_COORD(v1, 1);
        CLAMP_COORD(v2, 2);
        CLAMP_COORD(v3, 3);
        return get_nearest(v0, v1, v2, v3);
    }

    /* operator()() const performs a linear/bilinear/trilinear lookup from non-downsampled coordinates.
       Coordinates with downsampling factor of 1 are not interpolated. */
    
#define DO_LERP(d0, d1, frac) ((d0) + ((d1)-(d0))*(frac))
#define LERP_VARS(ndim) \
        int vi##ndim(v##ndim/downsample_dim##ndim); \
        int vip##ndim(MIN(vi##ndim+1, sizes_downsampled[ndim]-1)); \
        float vf##ndim((v##ndim%downsample_dim##ndim)/double(downsample_dim##ndim)); \
        ASSERT(vi##ndim >= 0 && vi##ndim < sizes_downsampled[ndim], "lookup out of bounds");

    INLINE real operator()(int v0) const {
        ASSERT(sizes.size() >= 1, "1D lookup in array of lesser dimensionality");
        if (downsample_dim0 == 1) {
            return get_nearest(v0);
        } else {
            LERP_VARS(0)
            
            real d0(data[vi0]);
            real d1(data[vip0]);
            return DO_LERP(d0, d1, vf0);
        }
    }

    INLINE real operator()(int v0, int v1) const {
        ASSERT(sizes.size() >= 2, "2D lookup in array of lesser dimensionality");
//        if (PRINT_COORDS(v0, v1)) { printf("array operator()\n"); }
        if (downsample_dim0 == 1 && downsample_dim1 == 1) {
//            if (PRINT_COORDS(v0, v1)) { printf("operator(), get nearest\n"); }
            return get_nearest(v0, v1);
        } else if (downsample_dim0 == 1) {
//            if (PRINT_COORDS(v0, v1)) { printf("operator(), lerp dim 1 (x)\n"); }
            LERP_VARS(1)
            int v0_stride = v0*stride[0];
            real d0(data[v0_stride+ vi1]);
            real d1(data[v0_stride+vip1]);
            return DO_LERP(d0, d1, vf1);
        } else if (downsample_dim1 == 1) {
//            if (PRINT_COORDS(v0, v1)) { printf("operator(), lerp dim 0 (y)\n"); }
            LERP_VARS(0)
            int v1_stride = v1*stride[1];
            real d0(data[vi0*stride[0] + v1_stride]);
            real d1(data[vip0*stride[0]+ v1_stride]);
            return DO_LERP(d0, d1, vf0);
        } else {
//            if (PRINT_COORDS(v0, v1)) { printf("operator(), lerp both dims\n"); }
            LERP_VARS(0)
            LERP_VARS(1)
            
            int vi0_stride = vi0 *stride[0];
            int vip0_stride = vip0*stride[0];
            
            real d00(data[vi0_stride + vi1]);
            real d01(data[vi0_stride +vip1]);
            real d10(data[vip0_stride+ vi1]);
            real d11(data[vip0_stride+vip1]);
            
            real d0(DO_LERP(d00, d01, vf1));
            real d1(DO_LERP(d10, d11, vf1));
            return DO_LERP(d0, d1, vf0);
        }
    }

#define LERP_2D_PREFIX(d0, d1) \
                LERP_VARS(d0) \
                LERP_VARS(d1) \
                \
                int vi##d0##_stride = vi##d0 *stride[d0]; \
                int vip##d0##_stride = vip##d0*stride[d0]; \
                int vi##d1##_stride = vi##d1*stride[d1]; \
                int vip##d1##_stride = vip##d1*stride[d1]; \

#define LERP_2D_SUFFIX(d0, d1) \
                real c00(datap[vi##d0##_stride + vi##d1##_stride ]); \
                real c01(datap[vi##d0##_stride + vip##d1##_stride]); \
                real c10(datap[vip##d0##_stride+ vi##d1##_stride ]); \
                real c11(datap[vip##d0##_stride+ vip##d1##_stride]); \
                \
                real c0(DO_LERP(c00, c01, vf##d1)); \
                real c1(DO_LERP(c10, c11, vf##d1)); \
                return DO_LERP(c0, c1, vf##d0);

    INLINE real operator()(int v0, int v1, int v2) const {
        ASSERT(sizes.size() >= 3, "3D lookup in array of lesser dimensionality");
       
        //printf("\nYES WE BE HERE!!!!!!!!!!!!!!!!!!!!!\n"); 
        if (downsample_dim2 == 1) {
            ASSERT(v2 >= 0 && v2 < sizes_downsampled[2], "3D lookup out of bounds on dimension 2");

            if (downsample_dim0 == 1 && downsample_dim1 == 1) {
                return get_nearest(v0, v1, v2);
            } else if (downsample_dim0 == 1) {
                LERP_VARS(1)
                int v0_stride = v0*stride[0];
                real d0(data[v0_stride+ vi1*stride[1]+v2]);
                real d1(data[v0_stride+vip1*stride[1]+v2]);
                return DO_LERP(d0, d1, vf1);
            } else if (downsample_dim1 == 1) {
                LERP_VARS(0)
                int v1_stride = v1*stride[1];
                real d0(data[vi0*stride[0] + v1_stride+v2]);
                real d1(data[vip0*stride[0]+ v1_stride+v2]);
                return DO_LERP(d0, d1, vf0);
            } else {
                LERP_VARS(0)
                LERP_VARS(1)
                
                int vi0_stride = vi0 *stride[0];
                int vip0_stride = vip0*stride[0];
                int vi1_stride = vi1*stride[1];
                int vip1_stride = vip1*stride[1];
                
                real d00(data[vi0_stride + vi1_stride +v2]);
                real d01(data[vi0_stride + vip1_stride+v2]);
                real d10(data[vip0_stride+ vi1_stride +v2]);
                real d11(data[vip0_stride+ vip1_stride+v2]);
                
                real d0(DO_LERP(d00, d01, vf1));
                real d1(DO_LERP(d10, d11, vf1));
                return DO_LERP(d0, d1, vf0);
            }
        } else {
            //printf("\n\n\n\n\n\nasdasdasdasdNOOOOOOOOOOOOOOOO\nwe hit the trilinear ASSERT\n\n\n");
            ASSERT2(0, "trilinear interpolation not implemented");
        }
    }

    INLINE real operator()(int v0, int v1, int v2, int v3) const {
        ASSERT(sizes.size() >= 4, "4D lookup in array of lesser dimensionality");
        
        if (downsample_dim0 == 1 && downsample_dim1 == 1 && downsample_dim2 == 1 && downsample_dim3 == 1) {
            return get_nearest(v0, v1, v2, v3);
        } else if (downsample_dim1 == 1 && downsample_dim2 == 1 && downsample_dim3 == 1) {             /* Lerp dim0 */
            LERP_VARS(0)
            real *datap = data + v1*stride[1] + v2*stride[2] + v3;
            real d0(datap[vi0*stride[0]]);
            real d1(datap[vip0*stride[0]]);
            return DO_LERP(d0, d1, vf0);
        } else if (downsample_dim0 == 1 && downsample_dim2 == 1 && downsample_dim3 == 1) {             /* Lerp dim1 */
            LERP_VARS(1)
            real *datap = data + v0*stride[0] + v2*stride[2] + v3;
            real d0(datap[vi1*stride[1]]);
            real d1(datap[vip1*stride[1]]);
            return DO_LERP(d0, d1, vf1);
        } else if (downsample_dim0 == 1 && downsample_dim1 == 1 && downsample_dim3 == 1) {             /* Lerp dim2 */
            LERP_VARS(2)
            real *datap = data + v0*stride[0] + v1*stride[1] + v3;
            real d0(datap[vi2*stride[2]]);
            real d1(datap[vip2*stride[2]]);
            return DO_LERP(d0, d1, vf2);
        } else if (downsample_dim2 == 1 && downsample_dim3 == 1) {                                     /* Lerp dims 0, 1 */
            LERP_2D_PREFIX(0, 1);
            real *datap = data + v2*stride[2] + v3;
            LERP_2D_SUFFIX(0, 1);
        } else if (downsample_dim0 == 1 && downsample_dim3 == 1) {                                     /* Lerp dims 1, 2 */
            LERP_2D_PREFIX(1, 2);
            real *datap = data + v0*stride[0] + v3;
            LERP_2D_SUFFIX(1, 2);
        } else if (downsample_dim0 == 1 && downsample_dim1 == 1) {                                     /* Lerp dims 2, 3 */
            LERP_2D_PREFIX(2, 3);
            real *datap = data + v0*stride[0] + v1*stride[1];
            LERP_2D_SUFFIX(2, 3);
        } else if (downsample_dim3 == 1) {                                                             /* Lerp dims 0, 1, 2 */
            LERP_VARS(0)
            LERP_VARS(1)
            LERP_VARS(2)
            
            int vi0_stride = vi0 *stride[0];
            int vip0_stride = vip0*stride[0];
            int vi1_stride = vi1*stride[1];
            int vip1_stride = vip1*stride[1];
            int vi2_stride = vi2*stride[2];
            int vip2_stride = vip2*stride[2];
            
            real *datap = data + v3;
            
            real d000(datap[vi0_stride + vi1_stride +vi2_stride]);
            real d010(datap[vi0_stride + vip1_stride+vi2_stride]);
            real d100(datap[vip0_stride+ vi1_stride +vi2_stride]);
            real d110(datap[vip0_stride+ vip1_stride+vi2_stride]);
            real d001(datap[vi0_stride + vi1_stride +vip2_stride]);
            real d011(datap[vi0_stride + vip1_stride+vip2_stride]);
            real d101(datap[vip0_stride+ vi1_stride +vip2_stride]);
            real d111(datap[vip0_stride+ vip1_stride+vip2_stride]);
            
            real d00(DO_LERP(d000, d010, vf1));
            real d10(DO_LERP(d100, d110, vf1));
            real d0 = DO_LERP(d00, d10, vf0);
            real d01(DO_LERP(d001, d011, vf1));
            real d11(DO_LERP(d101, d111, vf1));
            real d1 = DO_LERP(d01, d11, vf0);
            return DO_LERP(d0, d1, vf2);
        } else {
            ASSERT2(0, "4D linear interpolation requested is not implemented");
        }
    }
    
    INLINE real operator()(const vector<int> &index) const {
        return get_nearest(index);
    }

    
#if ASSIGN_ARRAY
    /* Assignment to array (A(y, x) = b) is disabled, use A.get_nearest(y, x) = b instead. */
    INLINE real &operator()(int v0) {
        return get_nearest(v0);
    }

    INLINE real &operator()(int v0, int v1) {
        return get_nearest(v0, v1);
    }

    INLINE real &operator()(int v0, int v1, int v2) {
        return get_nearest(v0, v1, v2);
    }

    INLINE real &operator()(int v0, int v1, int v2, int v3) {
        return get_nearest(v0, v1, v2, v3);
    }
    
    INLINE real &operator()(const vector<int> &index) {
        return get_nearest(index);
    }
    
#endif

    void clear() {
        memset((void *) data, 0, sizeof(real)*nelems);
    }
    
    void clear(real v) {
        for (int i = 0; i < nelems; i++) {
            data[i] = v;
        }
    }
    
    string str() const {
        string ans("[");
        if (dimensions() == 1) {
            for (int i = 0; i < nelems; i++) {
                ans += num_to_str(data[i]);
                if (i < nelems - 1) { ans += ", "; }
            }
        } else if (dimensions() == 2) {
            for (int y = 0; y < height(); y++) {
                ans += "[";
                for (int x = 0; x < width(); x++) {
                    ans += num_to_str((*this)(y,x));
                    if (x < width()-1) { ans += ", "; }
                }
                ans += "]";
                if (y < height()-1) { ans += ",\n "; }
            }
        }else if (dimensions() == 3) {
            for (int z = 0; z < channels(); z++) {
                ans += "\nChannel "+num_to_str(z)+":\n ";
                for (int y = 0; y < height(); y++) {
                    ans += "[";
                    for (int x = 0; x < width(); x++) {
                        ans += num_to_str((*this)(y,x,z));
                        if (x < width()-1) { ans += ", "; }
                    }
                    ans += "]";
                    if (y < height()-1) { ans += ",\n "; }
                }
                ans += "\n";
            }
        } else {
            ASSERT(false, "print not implement for dims > 3");
        }
        ans += "]";
        return ans;
    }

    static Array<real> random(const vector<int> sizes) {
        Array<real> ans(sizes);
        for (int i = 0; i < ans.nelems; i++) {
            ans.data[i] = rand_f();
        }
        return ans;
    }

//    CImg<real> getCImg(){ // Returns a grayscale CImg object
//    	ASSERT(sizes.size() >= 2, "need at least 2 dims to get image");
//		CImg<real> output(width(),height(),1,1,0);
//		for (int y = 0; y < height() ; y++){
//			for (int x = 0; x < width() ; x++){
//				output(x,y,0,0) = (*this)(y,x);
//			}
//		}
//		return output;
//    }

#if USE_HALIDE
    Halide::Image<real> getHalideImage(){
    	ASSERT(sizes.size() >= 2, "need at least 2 dims to get image");
		Halide::Image<real> output(width(),height(),3);
	    for (int y = 0; y < height() ; y++){
		    for (int x = 0; x < width() ; x++){
        		for (int z = 0; z < channels() ; z++){
			        output(x,y,z) = (*this)(y,x,z);
		        }
		    }
	    }
		return output;
    }
#endif

	real sum(){ // returns sum of all elements in the array
		real total=0;
        for (int i = 0; i < nelems; i++) {
            total += data[i];
        }
        return total;
    }
    
	void divideBy(real divisor){ // divides all elements in array by divisor
    	ASSERT(divisor != 0, "cannot divide by zero");
        #pragma omp parallel for
        for (int i = 0; i < nelems; i++) {
            data[i] /= divisor;
        }
    }
    
	void multiplyBy(real factor){ // multiplies all elements in array by factor
        #pragma omp parallel for
        for (int i = 0; i < nelems; i++) {
            data[i] *= factor;
        }
    }
    
	real getMax(){ 
		real currentMax = data[0];
        for (int i = 1; i < nelems; i++) {
        	currentMax = MAX(currentMax, data[i]);
        }
        return currentMax;
    }
    
	real getMin(){ 
		real currentMin = data[0];
        for (int i = 1; i < nelems; i++) {
        	currentMin = MIN(currentMin, data[i]);
        }
        return currentMin;
    }
    
	void rgb2gray(){ // Similar to rgb2gray in matlab
    	ASSERT(sizes.size() >= 3, "need at least 3 dims to convert to grayscale");
    	
    	real* data_grayscale = new real[height() * width()];
	    for (int y = 0; y < height() ; y++){
		    for (int x = 0; x < width() ; x++){
		        real cum_sum = 0;
		        for (int z = 0; z < channels(); z++){
                    real channel_weight = 1.0/channels();
                    if (channels() == 3) {
                        if (z == 0) { channel_weight = 0.299; }
                        else if (z == 1) { channel_weight = 0.5870; }
                        else if (z == 2) { channel_weight = 0.1140; }
                    }
				    cum_sum += get_nearest(y,x,z)*channel_weight;
			    }
			    data_grayscale[y*width()+x] = cum_sum;
		    }
	    }
        delete[] data;
        stride.resize(2);
        stride[0]=width();
        stride[1]=1;
        sizes.resize(2); // First two dimensions should be the same since we're just throwing away all the other dimensions
        data = data_grayscale;
    }

#define ARRAY_RGB2YUV() \
	    for (int y = 0; y < height() ; y+= downsample_dim0){ \
		    for (int x = 0; x < width() ; x+= downsample_dim1){ \
		        real old_R = get_nearest(y,x,R_COORD); \
		        real old_G = get_nearest(y,x,G_COORD); \
		        real old_B = get_nearest(y,x,B_COORD); \
                \
		        out.get_nearest(y,x,0) = 0.299   *old_R + 0.587   *old_G + 0.114   *old_B; \
		        out.get_nearest(y,x,1) = -0.14713*old_R + -0.28886*old_G + 0.436   *old_B; \
		        out.get_nearest(y,x,2) = 0.615   *old_R + -0.51499*old_G + -0.10001*old_B; \
		    } \
	    }

    void rgb2yuv(Array<real> &out) const {
        out.resize(this->sizes);

        #pragma omp parallel for
        ARRAY_RGB2YUV();
    }

    void rgb2yuv() {
        Array<real> &out(*this);

        #pragma omp parallel for
        ARRAY_RGB2YUV();
    }

    void yuv2rgb() {
        #pragma omp parallel for
	    for (int y = 0; y < height() ; y+= downsample_dim0){
		    for (int x = 0; x < width() ; x+= downsample_dim1){
		        real old_Y = get_nearest(y,x,0);
		        real old_U = get_nearest(y,x,1);
		        real old_V = get_nearest(y,x,2);

		        get_nearest(y,x,R_COORD) = old_Y + 1.13983*old_V;
		        get_nearest(y,x,G_COORD) = old_Y + -0.39465*old_U + -0.58060*old_V;
		        get_nearest(y,x,B_COORD) = old_Y + 2.03211*old_U;
		    }
	    }
    }
    
	void rgb2xyz(){ // code based on http://www.brucelindbloom.com/
#if DEBUG_TIME
        double T0 = wall_time();
#endif
    	ASSERT(sizes.size() == 3, "need at exactly 3 dims to convert from sRGB to XYZ");
    	ASSERT(getMax() <= 1.0, "need elements to be scaled, i.e. in the range [0.0, 1.0]");
    	ASSERT(getMin() >= 0.0, "need elements to be scaled, i.e. in the range [0.0, 1.0]");
    	
        static PowTable<real> pow_table(2.4);
        
        #pragma omp parallel for
        for (int i = 0; i < nelems; i++) {
//            data[i] = (data[i] <= 0.04045) ? data[i]/12.92 : pow(((data[i]+0.055)/1.055),2.4);
            data[i] = (data[i] <= 0.04045) ? data[i]/12.92 : pow_table(((data[i]+0.055)/1.055));
        }
    	double Tmid = wall_time();
        
        #pragma omp parallel for
	    for (int y = 0; y < height() ; y+= downsample_dim0){
		    for (int x = 0; x < width() ; x+= downsample_dim1){
		        real old_R = get_nearest(y,x,R_COORD);
		        real old_G = get_nearest(y,x,G_COORD);
		        real old_B = get_nearest(y,x,B_COORD);
		        /*
		        // D65 Coefficients
		        (*this)(y,x,X_COORD) = 0.4124564*old_R + 0.3575761*old_G + 0.1804375*old_B;
		        (*this)(y,x,Y_COORD) = 0.2126729*old_R + 0.7151522*old_G + 0.0721750*old_B;
		        (*this)(y,x,Z_COORD) = 0.0193339*old_R + 0.1191920*old_G + 0.9503041*old_B;
		        */
		        // D50 Coefficients
		        get_nearest(y,x,X_COORD) = 0.4360747*old_R + 0.3850649*old_G + 0.1430804*old_B;
		        get_nearest(y,x,Y_COORD) = 0.2225045*old_R + 0.7168786*old_G + 0.0606169*old_B;
		        get_nearest(y,x,Z_COORD) = 0.0139322*old_R + 0.0971045*old_G + 0.7141733*old_B;
		    }
	    }
#if DEBUG_TIME
        double T1 = wall_time();
        printf("rgb2xyz: %f seconds (%f %f)\n", T1-T0, Tmid-T0, T1-Tmid);
#endif
    }
    
	// We need a reference white for our conversions between XYZ and L*a*b*
	// We'll be using the default reference white that is used in the code snippet described here: http://www.mathworks.com/matlabcentral/answers/47601-how-to-convert-rgb-to-cie-lab-color-space
	// The snippet uses makecform: http://www.mathworks.com/help/images/ref/makecform.html
	// The default reference white used by makecform can be found here: http://www.mathworks.com/help/images/ref/whitepoint.html
	// This uses the ICC standard profile connection space illuminant, a 16-bit fractional approximation of D50
	
	#define REFERENCE_WHITE_X 0.9642;
	#define REFERENCE_WHITE_Y 1.0000;
	#define REFERENCE_WHITE_Z 0.8249;
	
	// As defined by the CIE standard
	#define CIE_EPSILON 0.008856
	#define CIE_KAPPA 903.3
	
	void xyz2lab(){ 
#if DEBUG_TIME
        double T0 = wall_time();
#endif
	    // code based on http://www.brucelindbloom.com/
    	ASSERT(sizes.size() == 3, "need at exactly 3 dims to convert from XYZ to L*a*b*");
    	
        static PowTable<real> pow_table(1.0/3.0);
        #pragma omp parallel for
	    for (int y = 0; y < height() ; y+= downsample_dim0){
		    for (int x = 0; x < width() ; x+= downsample_dim1){
            	real x_r = get_nearest(y,x,X_COORD) / REFERENCE_WHITE_X;
            	real y_r = get_nearest(y,x,Y_COORD) / REFERENCE_WHITE_Y;
            	real z_r = get_nearest(y,x,Z_COORD) / REFERENCE_WHITE_Z;
            	
//            	real f_x = (x_r > CIE_EPSILON) ? pow(x_r, 1.0/3.0) : (CIE_KAPPA*x_r+16)/116;
//            	real f_y = (y_r > CIE_EPSILON) ? pow(y_r, 1.0/3.0) : (CIE_KAPPA*y_r+16)/116;
//            	real f_z = (z_r > CIE_EPSILON) ? pow(z_r, 1.0/3.0) : (CIE_KAPPA*z_r+16)/116;
            	real f_x = (x_r > CIE_EPSILON) ? pow_table(x_r) : (CIE_KAPPA*x_r+16)/116;
            	real f_y = (y_r > CIE_EPSILON) ? pow_table(y_r) : (CIE_KAPPA*y_r+16)/116;
            	real f_z = (z_r > CIE_EPSILON) ? pow_table(z_r) : (CIE_KAPPA*z_r+16)/116;
            	
                get_nearest(y,x,L_COORD) = 116*f_y-16;
                get_nearest(y,x,A_COORD) = 500*(f_x-f_y);
                get_nearest(y,x,B_COORD) = 200*(f_y-f_z);
		    }
	    }
#if DEBUG_TIME
        double T1 = wall_time();
        printf("xyz2lab: %f seconds\n", T1-T0);
#endif
    }
    
	void rgb2lab(int scale_to_unity=0){
#if DEBUG_TIME
        double T0 = wall_time();
#endif
	    rgb2xyz();
	    xyz2lab();
        if (scale_to_unity) {
            multiplyBy(1.0/100.0);
        }
#if DEBUG_TIME
        double T1 = wall_time();
        printf("rgb2lab: %f seconds\n", T1-T0);
#endif
    }
    
	void lab2xyz(){
	    // code based on http://www.brucelindbloom.com/
    	ASSERT(sizes.size() == 3, "need at exactly 3 dims to convert from L*a*b* to XYZ");
    	
        #pragma omp parallel for
	    for (int y = 0; y < height() ; y+= downsample_dim0){
		    for (int x = 0; x < width() ; x+= downsample_dim1){
            	real f_y = (get_nearest(y,x,L_COORD)+16)/116;
            	real f_x = get_nearest(y,x,A_COORD)/500.0 + f_y;
            	real f_z = f_y - get_nearest(y,x,B_COORD)/200.0;
		        
            	real x_r = (f_x*f_x*f_x > CIE_EPSILON) ? f_x*f_x*f_x : (116*f_x-16)/CIE_KAPPA;
            	real y_r = (get_nearest(y,x,L_COORD) > CIE_KAPPA*CIE_EPSILON) ? f_y*f_y*f_y : get_nearest(y,x,L_COORD) / CIE_KAPPA;
            	real z_r = (f_z*f_z*f_z > CIE_EPSILON) ? f_z*f_z*f_z : (116*f_z-16)/CIE_KAPPA;
		        
                get_nearest(y,x,X_COORD) = x_r*REFERENCE_WHITE_X;
                get_nearest(y,x,Y_COORD) = y_r*REFERENCE_WHITE_Y;
                get_nearest(y,x,Z_COORD) = z_r*REFERENCE_WHITE_Z;
		    }
	    }
    }
    
	void xyz2rgb(){ // code based on http://www.brucelindbloom.com/
    	ASSERT(sizes.size() == 3, "need at exactly 3 dims to convert from XYZ to sRGB");
    	
        #pragma omp parallel for
	    for (int y = 0; y < height() ; y+= downsample_dim0){
		    for (int x = 0; x < width() ; x+= downsample_dim1){
		        real old_X = get_nearest(y,x,X_COORD);
		        real old_Y = get_nearest(y,x,Y_COORD);
		        real old_Z = get_nearest(y,x,Z_COORD);
		        
		        get_nearest(y,x,R_COORD) = 3.1338561*old_X + -1.6168667*old_Y + -0.4906146*old_Z;
		        get_nearest(y,x,G_COORD) = -0.9787684*old_X + 1.9161415*old_Y + 0.0334540*old_Z;
		        get_nearest(y,x,B_COORD) = 0.0719453*old_X + -0.2289914*old_Y + 1.4052427*old_Z;
		    }
	    }
	    
        static PowTable<real> pow_table(1.0/2.4);
        
        #pragma omp parallel for
        for (int i = 0; i < nelems; i++) {
///            data[i] = (data[i] <= 0.0031308) ? 12.92*data[i] : 1.055*pow(data[i],(1.0/2.4)) - 0.055;
            data[i] = (data[i] <= 0.0031308) ? 12.92*data[i] : 1.055*pow_table(data[i]) - 0.055;
        }
    	
    }
    
	void lab2rgb(int scale_to_unity=0){
#if DEBUG_TIME
        double T0 = wall_time();
#endif
        if (scale_to_unity) { multiplyBy(100.0); }
	    lab2xyz();
	    xyz2rgb();
#if DEBUG_TIME
        double T1 = wall_time();
        printf("lab2rgb: %f seconds\n", T1-T0);
#endif
	}
	
	void getChannel(int channel, Array<real> &result){
		result.resize({height(), width()});
	    for (int y = 0; y < height() ; y++){
		    for (int x = 0; x < width() ; x++){
		        result(y,x) = (*this)(y,x,channel);
		    }
	    }
    }
    
	void addThirdDimension(){
		// Changes the array (assuming it is 2D) into a 3D array of size height x width x 1
    	ASSERT(sizes.size() == 2, "need at exactly 2 dims to add third dimension");
        stride.resize(height(), width(), 1);
        nelems = 1;
        for (int i = stride.size()-1; i >= 0; i--) {
            stride[i] = nelems;
            nelems *= sizes[i];
        }
    }
	
	void clamp(real newMin=0.0, real newMax=1.0){ // Clamps (i.e. clips) values of array to being in the range [newMin,newMax]
        for (int i = 0; i < nelems; i++) {
            if (data[i] > newMax){
                data[i] = newMax;
            } else if (data[i] < newMin){
                data[i] = newMin;
            }
        }
    }
    
	void normalize(real newMin, real newMax){ // linearly remaps values of array from old range to new range [newMin,newMax]
		real oldMin = getMin();
		real oldMax = getMax();
        for (int i = 0; i < nelems; i++) {
	        //data[i] = (data[i]-oldMin) / (oldMax-oldMin) * (newMax-newMin) + newMin;
	        data[i] = newMax - (oldMax-data[i]) / (oldMax-oldMin) * (newMax-newMin);
        }
    }
    
    real max() {
		real ans = data[0];
        for (int i = 1; i < nelems; i++) {
	        ans = MAX(ans, data[i]);
        }
        return ans;
    }

    real min() {
		real ans = data[0];
        for (int i = 1; i < nelems; i++) {
	        ans = MIN(ans, data[i]);
        }
        return ans;
    }
    
	void normalize2(){ // Divides everything by the max, so now everything is less than 1.
        real oldMax(max());
        for (int i = 0; i < nelems; i++) {
	        data[i] =data[i] / oldMax;
        }
    }
    
    // Produces a sub-view of a given rectangle, may be outside the image's boundaries in which case accessing those pixels is undefined.
    const Array<real> *subview(int top, int left, int height, int width) const {
        ASSERT2(dimensions() == 3, "expected 3 dimension input image for subview()");
        Array<real> *ans = new Array<real>();
        delete ans->data;
        ans->dealloc = false;
        ans->data = &get_unsafe(top, left, 0);
        int nchannels = channels();
        ans->sizes = vector<int>({height, width, nchannels});
        ans->sizes_downsampled = ans->sizes;
        ans->nelems = height * width * nchannels;
        ans->stride = vector<int>({width*nchannels, nchannels, 1});
        
        return ans;
    }
    
	void toString(string& s){
    	ASSERT(sizes.size() <= 3, "toString not implemented for more than 3 dims");
	    s = "[";
	    int num_channels = (sizes.size() == 3) ? channels() : 1;
        for (int z = 0; z < num_channels; z++){
            if (sizes.size() == 3) {
                s += "\nChannel: "+num_to_str(z)+"\n ";
            }
            for (int y = 0; y < height(); y++){
                if (sizes.size() > 1){
                    s += "[";
	                for (int x = 0; x < width(); x++){
	                    if (sizes.size() > 2) {
	                        s += num_to_str(this->get_nearest(y,x,z));
                        } else {
		                    s += num_to_str(this->get_nearest(y,x));
                        }
	                    if (x != width()-1) { 
	                        s += ", ";
                        }
	                }
                    s += "]";
                    if (y != height()-1) { 
                        s += ",\n ";
                    }
                } else {
                    s += num_to_str(this->get_nearest(y));
                    if (y != height()-1) { 
                        s += ", ";
                    };
                }
            }
            if (sizes.size() == 3) {
                s += "\n";
            }
        }
        s += "]";
    }
    
#if ARRAY_USE_OPENCV
    template<class scalar>
    cv::Mat to_cv_vectorized() const {
        int type = boost::mpl::at<cvTypeChannels3, scalar>::type::value;
        ASSERT2(dimensions() == 3, "expected 3 dimensions for to_cv_vectorized");
        cv::Mat ans(2, &sizes[0], type, cv::Scalar::all(0));
        
        for (int y = 0; y < height(); y++) {
            scalar *row = ans.ptr<scalar>(y);
            for (int x = 0; x < width(); x++) {
                scalar *pixel = &row[x*3];
                real *this_pixel = data + y * stride[0] + x * stride[1];
                pixel[0] = (*this_pixel)[0];
                pixel[1] = (*this_pixel)[1];
                pixel[2] = (*this_pixel)[2];
            }
        }
        return ans;
    }
    
    /* Convert to OpenCV Mat */
    cv::Mat to_cv() const {
        int type = CV_TYPE(*this, real);
        if (type == UNKNOWN_CV_TYPE) { fprintf(stderr, "Cannot find OpenCV type for array, sizes=%s\n", vector_to_str_int<int>(sizes).c_str()); exit(1); }
        
        if (dimensions() != 2 && dimensions() != 3) {
            fprintf(stderr, "to_cv: dimensions not 2 or 3, cannot convert to cv::Mat\n"); exit(1);
        }
        cv::Mat ans(2, &sizes[0], type, cv::Scalar::all(0));
        
        if (dimensions() == 2) {
            for (int y = 0; y < height(); y++) {
                real *row = ans.ptr<real>(y);
                for (int x = 0; x < width(); x++) {
                    row[x] = (*this)(y, x);
                }
            }
        } else if (dimensions() == 3) {
            for (int y = 0; y < height(); y++) {
                real *row = ans.ptr<real>(y);
                for (int x = 0; x < width(); x++) {
                    real *pixel = &row[x*channels()];
                    for (int c = 0; c < channels(); c++) {
                        pixel[c] = (*this)(y, x, c);
                    }
                }
            }
        }
        return ans;
    }

    /* Convert to uint8_t OpenCV Mat */
    cv::Mat to_cv8() const {
        int type = CV_TYPE(*this, uint8_t);
        if (type == UNKNOWN_CV_TYPE) { fprintf(stderr, "Cannot find OpenCV type for array, sizes=%s\n", vector_to_str_int<int>(sizes).c_str()); exit(1); }
        
        if (dimensions() != 2 && dimensions() != 3) {
            fprintf(stderr, "to_cv: dimensions not 2 or 3, cannot convert to cv::Mat\n"); exit(1);
        }
        cv::Mat ans(2, &sizes[0], type, cv::Scalar::all(0));
        
        if (dimensions() == 2) {
            for (int y = 0; y < height(); y++) {
                uint8_t *row = ans.ptr<uint8_t>(y);
                for (int x = 0; x < width(); x++) {
                    row[x] = uint8_t((*this)(y, x)*255.9999);
                }
            }
        } else if (dimensions() == 3) {
            for (int y = 0; y < height(); y++) {
                uint8_t *row = ans.ptr<uint8_t>(y);
                for (int x = 0; x < width(); x++) {
                    uint8_t *pixel = &row[x*channels()];
                    for (int c = 0; c < channels(); c++) {
                        pixel[c] = uint8_t((*this)(y, x, c)*255.9999);
                    }
                }
            }
        }
        return ans;
    }
    
    Array(const cv::Mat &M, int dims=3) :data(NULL), dealloc(true)  {
        if (dims == 3) {
            resize(M.rows, M.cols, M.channels());
            for (int y = 0; y < height(); y++) {
                const real *row = M.ptr<real>(y);
                for (int x = 0; x < width(); x++) {
                    const real *pixel = &row[x*channels()];
                    for (int c = 0; c < channels(); c++) {
                        (*this)(y, x, c) = pixel[c];
                    }
                }
            }
        } else if (dims == 2) {
            resize(M.rows, M.cols);
            for (int y = 0; y < height(); y++) {
                const real *row = M.ptr<real>(y);
                for (int x = 0; x < width(); x++) {
                    (*this)(y, x) = row[x];
                }
            }
        } else {
            fprintf(stderr, "Array(cv::Mat, dims), dims is not 2 or 3, cannot convert from cv::Mat\n"); exit(1);
        }
    }
    
#endif
};

template<class real>
void save_color_image(const Array<real> &a, const char *filename);

#define MAX_LINE 65536

/* Save 2D array to ASCII file without any header. Compatible with MATLAB/numpy. */
template<class real>
void save_matrix(const Array<real> &A, string filename) {
    FILE *f = fopen(filename.c_str(), "wt");
    if (!f) { fprintf(stderr, "Error in save_matrix: Could not open %s\n", filename.c_str()); exit(1); }
    
    for (int y = 0; y < A.height(); y++) {
        for (int x = 0; x < A.width(); x++) {
            fprintf(f, "%.16g ", to_double(A(y, x)));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

/* Load 2D array from ASCII file without any header (loadAscii uses width height depth for header). Compatible with MATLAB/numpy. */
template<class real>
Array<real> load_matrix(string filename, bool verbose=false) {
    Array<real> ans;
    FILE *f = fopen(filename.c_str(), "rt");
    if (!f) { fprintf(stderr, "Error in load_matrix: Could not read %s\n", filename.c_str()); exit(1); }
    
    char buf[MAX_LINE];
    int rows = 0, cols = 0;
    while (1) {
        if (fgets(buf, MAX_LINE, f) == NULL) {
            break;
        }
        string s(strip(string(buf)));
        if (s.size() > 0) {
            rows++;
            string sp;
            while (1) {
                sp = str_replace(s, "  ", " ");
                if (sp == s) { break; }
                s = sp;
            }
            cols = count(s, ' ')+1;
        }
    }
    //exit(1);
    fseek(f, 0, SEEK_SET);
    
    ans.resize(rows, cols);
    
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            double value = 0.0;
            if (fscanf(f, "%lf", &value) != 1) {
                fprintf(stderr, "load_matrix: Error reading %s\n", filename.c_str());
                exit(1);
            }
            ans.get_nearest(row, col) = value;
        }
    }
    if (verbose) {
        printf("load_matrix: read %dx%d\n", rows, cols);
        printf("%s\n", ans.str().c_str());
    }
    fclose(f);
    return ans;
}

typedef boost::mpl::map<
  boost::mpl::pair<double,   boost::mpl::integral_c<int,      1> >,
  boost::mpl::pair<float,    boost::mpl::integral_c<int,      1> >,
  boost::mpl::pair<uint8_t,  boost::mpl::integral_c<uint8_t,  255> >,
  boost::mpl::pair<uint16_t, boost::mpl::integral_c<uint16_t, 65535> >,
  boost::mpl::pair<uint32_t, boost::mpl::integral_c<uint32_t, 65535> >,
  boost::mpl::pair<uint64_t, boost::mpl::integral_c<uint64_t, 65535> >,
  boost::mpl::pair<int8_t,   boost::mpl::integral_c<int8_t,   127> >,
  boost::mpl::pair<int16_t,  boost::mpl::integral_c<int16_t,  32767> >,
  boost::mpl::pair<int32_t,  boost::mpl::integral_c<int32_t,  65535> >,
  boost::mpl::pair<int64_t,  boost::mpl::integral_c<int64_t,  65535> > > imageMaxByType;

#define IMAGE_MAX_BY_TYPE(T) boost::mpl::at<imageMaxByType, T>::type::value
#define LOAD_SCALE_FACTOR(T) (IMAGE_MAX_BY_TYPE(T)/(255.0))

string get_additional_filename(const char *filename, int n);

bool string_ends_with(const char *a, const char *b);

template<class real>
Array<real> load_color_image(const char *filename, bool extra_channels=false, double scale=-1) {
//    CImg<real> image(filename);
    ASSERT(extra_channels==false, "Extra channels have not been implemented yet in the new png loader that doesn't depend on CImg");
    bool is_text = string_ends_with(filename, ".txt"); //strstr(filename, ".txt");
    if (is_text) {
        Array<real> ans;
        ans.loadAscii(string(filename));
        return ans;
    }
    
    bool is_pfm = string_ends_with(filename, ".pfm"); //strstr(filename, ".pfm");
    float *depth = NULL;
    int* w = new int;
    int* h = new int;
    unsigned depth_length;
    png::image< png::rgb_pixel > I;
    if (is_pfm) { 
        scale = 1.0;
        depth = read_pfm_file(filename, w, h);
        if (depth == NULL) {
            depth = read_pfm_file3(filename, w, h);
        }
        ASSERT(depth!=NULL, "Could not load .pfm");
        depth_length = (*w)*(*h)*3;
    } else {
        I = png::image< png::rgb_pixel >(filename);
        (*h) = I.get_height();
        (*w) = I.get_width();
    }
    Array<real> ans((*h), (*w), extra_channels ? 6: 3 /* we choose to ignore the alpha channel */ );
    if (scale < 0) { scale = LOAD_SCALE_FACTOR(real); }
    for (int y = 0; y < (int) ans.height(); y++) {
        int yflip = ans.height()-1-y;
        for (int x = 0; x < (int) ans.width(); x++) {
            if (is_pfm) { 
                ans.get_nearest(y, x, 0) = depth[yflip*(*w)*3+x*3+0];
                ans.get_nearest(y, x, 1) = depth[yflip*(*w)*3+x*3+1];
                ans.get_nearest(y, x, 2) = depth[yflip*(*w)*3+x*3+2];
            } else {
                png::rgb_pixel px = I.get_pixel(x,y);
                
                ans.get_nearest(y, x, 0) = (real)(px.red*scale);
                ans.get_nearest(y, x, 1) = (real)(px.green*scale);
                ans.get_nearest(y, x, 2) = (real)(px.blue*scale);
            }
        }
    }
    if (is_pfm) { delete[] depth; }
//    if (extra_channels) {
//        ASSERT2(false, "load_color_image with extra channels not implemented");
//        string filename_p = get_additional_filename(filename, 1);
//        printf("load_color_image: reading extra channels from %s\n", filename_p.c_str());
//        CImg<real> image2(filename_p.c_str());
//        for (int y = 0; y < (int) ans.height(); y++) {
//            for (int x = 0; x < (int) ans.width(); x++) {
//                for (int c = 0; c < (int) MIN(image2.spectrum(), 3); c++) {
////                    printf("reading channel %d (extra_channels)\n", c);
//                    ans.get_nearest(y, x, c+3) = image2(x, y, 0, c)*scale;
//                }
//            }
//        }
//    }
	delete w;
	delete h;
    return ans;
}

template<class real>
void image_to_uint32(const Array<real> &a, Array<uint32_t> &b) {
    b.resize(a.height(), a.width());
    if (a.dimensions() == 3 && a.channels() == 3) {
        for (int y = 0; y < a.height(); y++) {
            real     *a_row = a.data + a.stride[0] * y;
            uint32_t *b_row = b.data + b.stride[0] * y;
            for (int x = 0; x < a.width(); x++) {
                real     *a_pixel = a_row + a.stride[1] * x;
                uint32_t *b_pixel = b_row + b.stride[1] * x;
                int r = a_pixel[0] * 255;
                int g = a_pixel[1] * 255;
                int b = a_pixel[2] * 255;
                *b_pixel = (r) | ((g)<<8) | ((b) << 16);
            }
        }
    } else {
        fprintf(stderr, "error: expected a dimensions to be 3 and a channels to be 3\n"); ASSERT2(false, "bad number of dimensions or channels");
    }
}

template<class real>
void check_valid_image(const Array<real> &a, const char *name, int channels=-1) {
    if (channels < 0) { channels = a.channels(); }
    if (a.dimensions() == 3) {
        for (int y = 0; y < a.height(); y++) {
            for (int x = 0; x < a.width(); x++) {
                for (int c = 0; c < channels; c++) {
                    if (std::isinf(a.get_nearest(y, x, c))) {
                        fprintf(stderr, "inf value at %d, %d, %d in %s\n", x, y, c, name); exit(1);
                    }
                    if (std::isnan(a.get_nearest(y, x, c))) {
                        fprintf(stderr, "nan value at %d, %d, %d in %s\n", x, y, c, name); exit(1);
                    }
                }
            }
        }
    }
}

template<class real>
void save_color_image(const Array<real> &a, const char *filename) {
//    printf("In save_color_image\n"); fflush(stdout);
    if (string_ends_with(filename, ".txt")) {
        a.saveAscii(string(filename));
        return;
    }
    
    if (a.dimensions() == 4) {
        Array<real> slice;
        slice.resize(a.sizes[1], a.sizes[2], a.sizes[3]);
        for (int z = 0; z < a.sizes[0]; z++) {
            for (int y = 0; y < a.sizes[1]; y++) {
                for (int x = 0; x < a.sizes[2]; x++) {
                    for (int c = 0; c < a.sizes[3]; c++) {
                        slice.get_nearest(y, x, c) = a.get_nearest(z, y, x, c);
                    }
                }
            }
            string filename_prefix = filename;
            if (strstr(filename, ".")) {
                filename_prefix = filename_prefix.substr(0, strstr(filename, ".")-filename);
            }
            char slice_filename[256];
            sprintf(slice_filename, "%s_%02d.pfm", filename_prefix.c_str(), z);
            save_color_image(slice, slice_filename);
        }
        return;
    }
    ASSERT2(a.dimensions() == 3, "expected 3D image in save_color_image");
    
//    printf("save_color_image, making png::image\n"); fflush(stdout);
//    CImg<real> image(a.width(), a.height(), 1, MIN(a.channels(), 3));
    png::image< png::rgb_pixel > I(a.width(), a.height());
    real max_value(255);
    bool is_pfm = string_ends_with(filename, ".pfm");
    double scale = is_pfm ? 1.0: 1.0/LOAD_SCALE_FACTOR(real);
    int num_channels = MIN(a.channels(), 3);
    unsigned depth_length = a.width()*a.height()*num_channels;
    float *depth = is_pfm ? new float[depth_length] : NULL;
    float *colors = new float[num_channels];
    
//    printf("save_color_image, created depth and colors, writing data to array\n"); fflush(stdout);
    for (int y = 0; y < (int) a.height(); y++) {
        for (int x = 0; x < (int) a.width(); x++) {
            for (int c = 0; c < num_channels; c++) {
                float value(a(y, x, c)*scale);
                if (!is_pfm) {
                    if (value < 0) { value = 0; }
                    else if (value > max_value) { value = max_value; }
                }
//                image(x, y, 0, c) = value;
                colors[c] = value;
            }
            ASSERT(num_channels==1 || num_channels==3, "save_color_image not implemented for images where the number of channels is not 1 or 3");
            if (!is_pfm) { // save png
                if (num_channels == 3) {
                    I[y][x] = png::rgb_pixel(CLAMP_INT(colors[0]), CLAMP_INT(colors[1]), CLAMP_INT(colors[2]));
                } else if (num_channels == 1) {
                    I[y][x] = png::rgb_pixel(CLAMP_INT(colors[0]), CLAMP_INT(colors[0]), CLAMP_INT(colors[0]));
                }
            } else {
                int yflip = (a.height()-1-y);
                if (num_channels == 3) {
                    depth[a.width()*yflip*3+x*3+0] = (float)(colors[0]);
                    depth[a.width()*yflip*3+x*3+1] = (float)(colors[1]);
                    depth[a.width()*yflip*3+x*3+2] = (float)(colors[2]);
                } else if (num_channels == 1) {
                    depth[a.width()*yflip+x] = (float)(colors[0]);
                }
            }
        }
    }
//    printf("save_color_image, before actual save\n"); fflush(stdout);
    if (is_pfm) { // save pfm
        if (num_channels == 3) {
            write_pfm_file3(filename, depth, a.width(), a.height());
        } else if (num_channels == 1) {
            write_pfm_file(filename, depth, a.width(), a.height());
        }
    } else {
//    image.save(filename);
        I.write(filename);
    }
    
//    printf("save_color_image, checking for additional channels\n"); fflush(stdout);
    if (is_pfm && a.channels() > 3) { depth_length = a.width()*a.height()*3; }
    for (int count = 3;; count += 3) {
//        printf("save_color_image, saving additional channel if in pfm mode\n"); fflush(stdout);
        if (is_pfm && a.channels() > count) {
            string filename_p = get_additional_filename(filename, count/3+1);
            printf("save_color_image: saving to %s and %s (more than 3 channels)\n", filename, filename_p.c_str());
            for (int i = 0; i < depth_length; i++) { depth[i] = 0; }
            for (int y = 0; y < (int) a.height(); y++) {
                int yflip = (a.height()-1-y);
                for (int x = 0; x < (int) a.width(); x++) {
                    for (int c = count; c < (int) MIN(a.channels(), count+3); c++) {
                        float value(a(y, x, c)*scale);
                        if (!is_pfm) {
                            if (value < 0) { value = 0; }
                            else if (value > max_value) { value = max_value; }
                        }
                        //image(x, y, 0, c-count) = value;
                        depth[a.width()*yflip*3+x*3+c-count] = (float)(value);
                    }
                }
            }
            //image.save(filename_p.c_str());
            write_pfm_file3(filename_p.c_str(), depth, a.width(), a.height());
        } else {
            break;
        }
    }
//    printf("save_color_image, deallocating\n"); fflush(stdout);
    delete[] colors;
//    printf("save_color_image, after delete colors\n"); fflush(stdout);
    delete[] depth;
//    printf("save_color_image, end of function\n"); fflush(stdout);
}

template<class real>
Array<real> zeros(const vector<int> sizes) {
    Array<real> ans(sizes);
    ans.clear();
    return ans;
}

template<class real>
Array<real> resample_nearest(const Array<real> &a, const vector<int> &sizes, bool isolate_center=false) {
    ASSERT(sizes.size() == 1 || sizes.size() == 2, "resample_nearest expected 1D or 2D input");
    if (sizes == a.sizes) { return a; }
    Array<real> ans(sizes);
    if (sizes.size() == 2) {
        int hp = ans.height(), wp = ans.width();
        int h = a.height(), w = a.width();
        for (int yp = 0; yp < hp; yp++) {
            int y = hp > 1 ? (yp*(h)/(hp)): 0;
            if (isolate_center) {
                if (yp == hp/2)   { y = h/2; }
                if (y == h/2) {
                    if (yp > hp/2) { y++; } else if (yp < hp/2) { y--; }
                }
/*
 if (yp == hp/2) { y = h/2; }
 else if (yp < hp/2) { y = hp/2 > 1 ? (yp*(h/2-1)/(hp/2-1)): 0; }
 else { y = hp/2 > 1 ? ((yp-hp/2-1)*(h/2-1)/(hp/2-1))+(h/2): h-1; }
 ASSERT(y >= 0 && y < h, "y out of bounds in resample_nearest");
*/
            }
 
            for (int xp = 0; xp < wp; xp++) {
                int x = wp > 1 ? (xp*(w)/(wp)): 0;
                if (isolate_center) {
                    if      (xp == wp/2)   { x = w/2; }
//                    else if (xp == wp/2-1) { x = w/2-1; }
//                    else if (xp == wp/2+1) { x = w/2+1; }
                    if (x == w/2) {
                        if (xp > wp/2) { x++; } else if (xp < wp/2) { x--; }
                    }
                }
                ans(yp, xp) = a(y, x);
            }
        }
    } else if (sizes.size() == 1) {
        int sp = ans.size();
        int s = a.size();
        for (int ip = 0; ip < sp; ip++) {
            int i = sp > 1 ? (ip*(s-1)/(sp-1)): 0;
            ans(i) = a(i);
        }
    } else {
        ASSERT(0, "resample_nearest expected 1D or 2D input");
    }
    return ans;
}

template<class real>
bool equals(const Array<real> &A, const Array<real> &B, real epsilon=0) {
    if (A.sizes != B.sizes) { return false; }
    for (int i = 0; i < A.nelems; i++) {
        real d(A.data[i]-B.data[i]);
        if (d < 0) { d = -d; }
        if (d > epsilon) { return false; }
    }
    return true;
}

/* filter_iir: Templated recursive (IIR) filter. Most of the implementation is in approx.py -- these are just some helper macros.
   Filters sequentially in directions +x, +y, -x, -y. Note +x/+y, -x/-y filtering occurs simultaneously in pairs.
   Filter is out[i] = a0*in[i]    + a1*in[i-1]  + a2*in[i-2]  + a3*in[i-3] +
                      b1*out[i-1] + b2*out[i-2] + b3*out[i-3] + b4*out[i-4].
   Zero coefficients are omitted from the computation. Modifies values in 'in'. */

#define make_zero0_1(real) 0
#define make_zero0_2(real) 0
#define make_zero0_3(real) 0
#define make_zero0_4(real) 0
#define make_zero2_1(real) { 0.0 }
#define make_zero2_2(real) { 0.0, 0.0 }
#define make_zero2_3(real) { 0.0, 0.0, 0.0 }
#define make_zero2_4(real) { 0.0, 0.0, 0.0, 0.0 }

#define IIR_IN_BOUNDS(xdir, ydir, delta) (ydir == 0 || (unsigned) (y-delta*ydir) < uh) && (xdir == 0 || (unsigned) (x-delta*xdir) < uw)

vector<float> gaussian_kernel(float sigma);

/* Crop or zero pad 2D input array A to new size w x h. */
template<class real>
Array<real> zero_pad_center(const Array<real> &A, int w, int h) {
    vector<int> sizes(w, h);
    Array<real> ans(sizes);
    int w0 = A.width();
    int h0 = A.height();
    ASSERT2(A.dimensions() == 2, "expected 2D array in zero_pad_center");
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int x0 = x-w/2+w0/2;
            int y0 = y-h/2+h0/2;
            if ((unsigned) x0 < (unsigned) w0 && (unsigned) y0 < (unsigned) h0) {
                ans.get_nearest(y, x) = A(y0, x0);
            } else {
                ans.get_nearest(y, x) = 0.0;
            }
        }
    }
    return ans;
}

/* Shift array by (dx, dy) with zero padding */
template<class real>
void shift_array(const Array<real> &A, Array<real> &shifted, int dx, int dy) {
    int h = A.height();
    int w = A.width();
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int xp = x-dx;
            int yp = y-dy;
            if ((unsigned) xp < (unsigned) w && (unsigned) yp < (unsigned) h) {
                shifted(y,x) = A(yp,xp);
            } else {
                shifted(y,x) = 0.0;
            }
        }
    }
}

#define DECLARE_OPERATORS(OPER) \
template<class real> \
Array<real> operator OPER (const Array<real> &a, const Array<real> &b) { \
    ASSERT(a.sizes == b.sizes, "a and b sizes differ for a * b"); \
    Array<real> ans(a.sizes); \
    for (int i = 0; i < a.nelems; i++) { \
        ans.data[i] = a.data[i] OPER b.data[i]; \
    } \
    return ans; \
} \
\
template<class real> \
Array<real> operator OPER (const Array<real> &a, real b) { \
    Array<real> ans(a.sizes); \
    for (int i = 0; i < a.nelems; i++) { \
        ans.data[i] = a.data[i] OPER b; \
    } \
    return ans; \
}

DECLARE_OPERATORS(*);

#define COPY_ARRAY(source) \
void copy_array(const Array<source> &src, Array<source> &dest);

COPY_ARRAY(double);
COPY_ARRAY(float);
COPY_ARRAY(uint8_t);
COPY_ARRAY(uint16_t);
COPY_ARRAY(uint32_t);
COPY_ARRAY(int8_t);
COPY_ARRAY(int16_t);
COPY_ARRAY(int32_t);

#undef COPY_ARRAY

template<class ArrayS, class ArrayT>
ArrayT cast_array(const ArrayS &source) {
    ArrayT dest;
    copy_array(source, dest);
    return dest;
}

#define XY_TO_INT(x, y) ((x)|((y)<<16))
#define INT_TO_X(v) ((v)&(65535))
#define INT_TO_Y(v) ((v)>>16)

#define COPY_FROM(I, xsrc, ysrc) \
    for (int channel = 0; channel < I.channels(); channel++) { \
        I.get_nearest(y, x, channel) = I.get_nearest(ysrc, xsrc, channel); \
    }

#define ADD_FRONTIER_NEIGHBORS() \
                if (x > 0   && !sample.get_nearest(y,x-1,0)) { frontier.insert(XY_TO_INT(x-1, y)); } \
                if (y > 0   && !sample.get_nearest(y-1,x,0)) { frontier.insert(XY_TO_INT(x, y-1)); } \
                if (x+1 < w && !sample.get_nearest(y,x+1,0)) { frontier.insert(XY_TO_INT(x+1, y)); } \
                if (y+1 < h && !sample.get_nearest(y+1,x,0)) { frontier.insert(XY_TO_INT(x, y+1)); }

template<class Array>
void copy_nearest_brute(Array &I, Array &sample) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif
    ASSERT2(I.dimensions() == 3, "I dimensions not 3");
    ASSERT2(sample.dimensions() == 3, "sample dimensions not 3");
    //ASSERT2(I.sizes == sample.sizes, "expected I.sizes to be sample.sizes");
    
    unordered_set<int> frontier;              /* Locations not sampled but with sampled neighbor. */
    int h = MIN(sample.height(), I.height());
    int w = MIN(sample.width(), I.width());
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (sample.get_nearest(y,x,0)) {
                ADD_FRONTIER_NEIGHBORS();
            }
        }
    }
    while (frontier.size()) {
        for (auto pos: frontier) {
            int x = INT_TO_X(pos);
            int y = INT_TO_Y(pos);
            ASSERT2(!sample.get_nearest(y,x,0), "sampled location should not be in frontier");
            if      (x > 0   && sample.get_nearest(y,x-1,0)) { COPY_FROM(I, x-1,y); }
            else if (y > 0   && sample.get_nearest(y-1,x,0)) { COPY_FROM(I, x,y-1); }
            else if (x+1 < w && sample.get_nearest(y,x+1,0)) { COPY_FROM(I, x+1,y); }
            else if (y+1 < h && sample.get_nearest(y+1,x,0)) { COPY_FROM(I, x,y+1); }
            else { ASSERT2(false, "no pixel to copy from in frontier"); }
            sample.get_nearest(y,x,0) = 1;
        }

        auto frontier_copy(frontier);
        frontier.clear();
        for (auto pos: frontier_copy) {
            int x = INT_TO_X(pos);
            int y = INT_TO_Y(pos);
            ADD_FRONTIER_NEIGHBORS();
        }
    }
#if DEBUG_TIME
    double T1 = wall_time();
    printf("copy_nearest_brute in %f seconds\n", T1-T0);
#endif
}

#define VBLUR_AT(x0, y0, c0) \
((current->get_clamp(y0-2, x0, c0)+current->get_clamp(y0+2, x0, c0))*coeff[0] + (current->get_clamp(y0-1, x0, c0)+current->get_clamp(y0+1, x0, c0))*coeff[1] + current->get_clamp(y0, x0, c0)*coeff[2])

#define GAUSSIAN_PYRAMID_COEFFS() \
    float coeff[] = { 1.0/16.0, 4.0/16.0, 6.0/16.0, 4.0/16.0, 1.0/16.0 }; \
    const int Kw = 5;

#define UPDATE_VBLUR(x, dx) \
                for (int c = 0; c < channels; c++) { \
                    vblur[(dx)+2][c] = VBLUR_AT(((x)+(dx))*2, y*2, c); \
                }

template<class Array, class real, class atom, int channels>
vector<Array *> create_gaussian_pyramid(Array &input) {
    static map<Array *, vector<void *> > reusable_mem;
    vector<void *> reusable_vec;
    bool has_reusable = false;
    if (reusable_mem.count(&input)) {
        reusable_vec = reusable_mem[&input];
        has_reusable = true;
    }
    
    Array *current = &input;
    vector<Array *> pyr({current});
    
    GAUSSIAN_PYRAMID_COEFFS();
    
    for (int i = 0;; i++) {
        vector<int> sizes_next({current->height()/2, current->width()/2, channels});
        if (sizes_next[0] < 1 || sizes_next[1] < 1) { break; }
        
        Array *next;
        if (has_reusable) {
            ASSERT2(i < (int) reusable_vec.size(), "i out of bounds in reusable memory");
            next = new Array(sizes_next, reusable_vec[i]);
        } else {
            next = new Array(sizes_next);
            reusable_vec.push_back((void *) next->data);
            next->dealloc = false;
        }
        
        #pragma omp parallel for
        for (int y = 0; y < next->height(); y++) {
            atom p[channels*Kw];
            atom *vblur[5] = { &p[0], &p[channels], &p[channels*2], &p[channels*3], &p[channels*4] };

            for (int dx = -2; dx <= 2; dx++) {
                UPDATE_VBLUR(0, dx);
            }

            for (int x = 0; x < next->width(); x++) {
                UPDATE_VBLUR(x, 2);

                for (int c = 0; c < channels; c++) {
                    next->get_nearest(y, x, c) = (vblur[0][c] + vblur[4][c]) * coeff[0] + (vblur[1][c] + vblur[3][c]) * coeff[1] + vblur[2][c] * coeff[2];
                }

                for (int i = 0; i < Kw-1; i++) {
                    for (int c = 0; c < channels; c++) {
                        vblur[i][c] = vblur[i+1][c];
                    }
                }
            }
        }
        pyr.push_back(next);
        current = next;
    }
    
    if (!has_reusable) {
        reusable_mem[&input] = reusable_vec;
    }

    return pyr;
}

template<class Array, class real, class atom, int channels>
void save_gaussian_pyramid(const vector<Array *> &pyr) {
    for (int i = 0; i < (int) pyr.size(); i++) {
        char buf[256];
        sprintf(buf, "pyramid%02d.png", i);
        save_color_image(*pyr[i], buf);
    }
}

template<class Array, class real, class atom, int channels>
void free_gaussian_pyramid(const vector<Array *> &pyr) {
    for (int i = 1; i < (int) pyr.size(); i++) {
        delete pyr[i];
    }
}

template<class Array, class real, class atom, int channels>
void upsample_gaussian_pyramid(const Array &input, Array &output, double mul_existing=0.0, double mul_new=1.0) {
    GAUSSIAN_PYRAMID_COEFFS();
    for (int i = 0; i < 5; i++) { coeff[i] *= 2; }
    
#define UPDATE_VBLUR_EVEN(dx) \
for (int c = 0; c < channels; c++) { \
    vblur[dx+1][c] = (input.get_clamp(y2-1, x2+dx, c) + input.get_clamp(y2, x2+dx, c)) * coeff[0] + input.get_clamp(y2+1, x2+dx, c) * coeff[2]; \
}

#define UPDATE_VBLUR_ODD(dx) \
for (int c = 0; c < channels; c++) { \
    vblur[dx+1][c] = (input.get_clamp(y2, x2+dx, c) + input.get_clamp(y2+1, x2+dx, c)) * coeff[1]; \
}

#define UPDATE_VBLUR_ALL() \
            if (y & 1 == 0) { \
                for (int dx = -1; dx <= 1; dx++) { \
                    UPDATE_VBLUR_EVEN(dx); \
                } \
            } else { \
                for (int dx = -1; dx <= 1; dx++) { \
                    UPDATE_VBLUR_ODD(dx); \
                } \
            }

#define UPDATE_XLOOP() \
        for (int x = 0; x < output.width(); x += 2) { \
            int x2 = x/2; \
\
            if (y & 1 == 0) {         \
                UPDATE_VBLUR_EVEN(1); \
            } else {                  \
                UPDATE_VBLUR_ODD(1);  \
            } \
\
            for (int c = 0; c < channels; c++) {

#define UPDATE_END_XLOOP() \
            } \
\
            for (int i = 0; i < 2; i++) { \
                for (int c = 0; c < channels; c++) { \
                    vblur[i][c] = vblur[i+1][c]; \
                } \
            } \
        }

    /* Blur x and y simultaneously */
    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        atom p[channels*3];
        atom *vblur[3] = { &p[0], &p[channels], &p[channels*2] };
        int y2 = y/2;
        
        {
            int x2 = 0;
            UPDATE_VBLUR_ALL();
        }

        if (mul_existing == 0.0 && mul_new == 1.0) {
            UPDATE_XLOOP();
            output.get_nearest(y, x,   c) = (vblur[0][c] + vblur[2][c]) * coeff[0] + vblur[1][c] * coeff[2];    /* Even x */
            output.get_nearest(y, x+1, c) = (vblur[1][c] + vblur[2][c]) * coeff[1];                             /* Odd x */
            UPDATE_END_XLOOP();
        } else {
            UPDATE_XLOOP();
            output.get_nearest(y, x,   c) = output.get_nearest(y, x,   c) * mul_existing + ((vblur[0][c] + vblur[2][c]) * coeff[0] + vblur[1][c] * coeff[2]) * mul_new;    /* Even x */
            output.get_nearest(y, x+1, c) = output.get_nearest(y, x+1, c) * mul_existing + ((vblur[1][c] + vblur[2][c]) * coeff[1]) * mul_new;                             /* Odd x */
            UPDATE_END_XLOOP();
        }
    }
}

template<class Array, class real, class atom, int channels>
void gaussian_blur_gaussian_pyramid(Array &input, Array &output, double sigma) {
//    printf("gaussian_blur_gaussian_pyramid\n");
    ASSERT2(input.dimensions() == 3, "expected 3 dimensional input to gaussian_blur_gaussian_pyr");
    output.resize(input.sizes);
    double T0 = wall_time();
    auto pyr = create_gaussian_pyramid<Array, real, atom, channels>(input);
    double T1 = wall_time();
        
    if (sigma < 0.001) { sigma = 0.001; }
    real level = log2(sigma/1.06988594);//+1;
    int ilevel(level);
    real flevel = level-ilevel;
//    printf("original ilevel=%d, pyr.size=%d\n", ilevel, pyr.size());
    
    if (ilevel >= int(pyr.size()-1)) {
        ilevel = pyr.size()-2;
        flevel = 1.0;
    }
    if (ilevel < 0) {
        ilevel = 0;
        flevel = 0.0;
    }
//    printf("ilevel=%d, flevel=%f, level=%f\n", ilevel, flevel, level);
    
    real flevel_compl = 1-flevel;
    
//    output.resize(pyr[ilevel]->sizes);
//    printf("gaussian_blur_gaussian_pyramid (before first resample), level=%f, ilevel=%d, flevel=%f\n", level, ilevel, flevel);
//    upsample_gaussian_pyramid<Array, real, atom, channels>(*pyr[ilevel+1], output);
    pyr[0] = &output;
    if (ilevel == 0) {
        pyr[ilevel]->assign(input);
    }
    upsample_gaussian_pyramid<Array, real, atom, channels>(*pyr[ilevel+1], *pyr[ilevel], flevel_compl, flevel);
//    ASSERT2(current.sizes == pyr[ilevel].sizes, "sizes differs from Gaussian pyramid size");
/*    const Array *pyr_ilevel = pyr[ilevel];
    for (int y = 0; y < output.height(); y++) {
        for (int x = 0; x < output.width(); x++) {
            for (int c = 0; c < channels; c++) {
                output.get_nearest(y, x, c) = output.get_nearest(y, x, c)*flevel + pyr_ilevel->get_nearest(y, x, c)*flevel_compl;
            }
        }
    }*/
//    printf("gaussian_blur_gaussian_pyramid (done first resample)\n"); fflush(stdout);
    //output.assign(*pyr[ilevel]);
//    pyr[0] = &output;
    
    for (int i = ilevel; i > 0; i--) {
//        printf("gaussian_blur_gaussian_pyramid (making output_copy), i=%d\n", i); fflush(stdout);
//        Array output_copy;
//        output_copy.assign(output);
//        printf("gaussian_blur_gaussian_pyramid (resizing output)\n"); fflush(stdout);
//        output.resize(pyr[i-1]->sizes);
//        printf("gaussian_blur_gaussian_pyramid (upsampling)\n"); fflush(stdout);
        upsample_gaussian_pyramid<Array, real, atom, channels>(*pyr[i], *pyr[i-1]);
//        ASSERT2(output.sizes == pyr[i-1]->sizes, "sizes differs from Pyramid level i-1 size");
//        printf("gaussian_blur_gaussian_pyramid (done upsampling), i=%d\n", i); fflush(stdout);
    }

//    save_gaussian_pyramid<Array, real, atom, channels>(pyr);

//    printf("gaussian_blur_gaussian_pyramid (done computation)\n"); fflush(stdout);
    free_gaussian_pyramid<Array, real, atom, channels>(pyr);
    double T2 = wall_time();
    printf("%f %f\n", T1-T0, T2-T1);
//    printf("gaussian_blur_gaussian_pyramid (done free)\n"); fflush(stdout);
}

template<class Array, class real, class atom, int channels, class Array_out, int order>
void cumsum_xy(const Array &input, Array_out &output) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif
    output.resize(input.sizes);
    ASSERT2(output.dimensions() == 3, "expected 3 dimensions");
    
    /* Cumulative sum x, repeated 'order' times. */
    
    /* Difference equations for order 2:
        x(i)                                          (Input)

        y(i) = x(i) + y(i-1),      y(1) = x(1)        (Cumulative sum of order 1)

        z(i) = y(i) + z(i-1),      z(1) = y(1)        (Cumulative sum of order 2)
                    Thus,
                    y(i) = z(i) - z(i-1)

        z(i) = x(i) + y(i-1) + z(i-1)
              = x(i) + z(i-1) - z(i-2) + z(i-1)
              = x(i) + 2*z(i-1) - z(i-2)

        z(2) = y(2) + z(1)
                = x(2) + y(1) + z(1)
                = x(2) * 2 * z(1) */

    
    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        for (int c = 0; c < channels; c++) {
            output.get_nearest(y, 0, c) = input.get_nearest(y, 0, c);
        }
        if (order == 1) {
            for (int x = 1; x < output.width(); x++) {
                for (int c = 0; c < channels; c++) {
                    output.get_nearest(y, x, c) = output.get_nearest(y, x-1, c) + input.get_nearest(y, x, c);
                }
            }
        } else if (order == 2) {
            for (int c = 0; c < channels; c++) {
                output.get_nearest(y, 1, c) = input.get_nearest(y, 1, c) + 2 * output.get_nearest(y, 0, c);
            }
            for (int x = 2; x < output.width(); x++) {
                for (int c = 0; c < channels; c++) {
                    output.get_nearest(y, x, c) = 2*output.get_nearest(y, x-1, c) + input.get_nearest(y, x, c) - output.get_nearest(y, x-2, c);
                }
            }
        } else {
            ASSERT2(false, "order > 2 not implemented");
        }
    }
#if DEBUG_TIME
    double Tmid = wall_time();
#endif

    /* Cumulative sum y */
    #pragma omp parallel for
    for (int x = 0; x < output.width(); x++) {
        if (order == 1) {
            for (int y = 1; y < output.height(); y++) {
                for (int c = 0; c < channels; c++) {
                    output.get_nearest(y, x, c) += output.get_nearest(y-1, x, c);
                }
            }
        } else if (order == 2) {
            for (int c = 0; c < channels; c++) {
                output.get_nearest(1, x, c) += 2 * output.get_nearest(0, x, c);
            }

            for (int y = 2; y < output.height(); y++) {
                for (int c = 0; c < channels; c++) {
                    output.get_nearest(y, x, c) += 2*output.get_nearest(y-1, x, c) - output.get_nearest(y-2, x, c);
                }
            }
        } else {
            ASSERT2(false, "order > 2 not implemented");
        }
    }

#if DEBUG_TIME
    double T1 = wall_time();
    printf("cumsum_xy: %f (%f x)\n", T1-T0, Tmid-T0);
#endif
}

#define OUTPUT_ORDER2(add13, add23, add33, add31, add32) \
                for (int c = 0; c < channels; c++) { \
                    output.get_nearest(y, x, c) = (   temp.get_nearest(y1,x1,c) - 2*temp.get_nearest(y1,x2,c) +   (temp.get_nearest(y1,x3,c)+(add13)) \
                                                   -2*temp.get_nearest(y2,x1,c) + 4*temp.get_nearest(y2,x2,c) - 2*(temp.get_nearest(y2,x3,c)+(add23)) \
                                                    + (temp.get_nearest(y3,x1,c)+(add31)) - 2*(temp.get_nearest(y3,x2,c)+(add32)) +   (temp.get_nearest(y3,x3,c)+(add33))) * scale; \
                }

#define GET_SPACING_SCALE() \
    double sigma = sigma_arr(y, x, 0) * 2; \
    int spacing = int(sigma); \
    if (spacing <= 1) { \
        for (int c = 0; c < channels; c++) { \
            output.get_nearest(y, x, c) = input.get_nearest(y, x, c); \
        } \
        continue; \
    } else if (spacing < 2) { spacing = 2; } \
    double scale; \
    if (order == 2) { \
        int spacing2 = spacing*spacing; \
        scale = 1.0/(spacing2*spacing2); \
    } else if (order == 1) { \
        scale = 1.0/(2*spacing*2*spacing); \
    }

/* From Heckbert 1986, Filtering by Repeated Integration. */
template<class ArrayT, class real, class atom, int channels, int order=2>
void gaussian_blur_sat(const ArrayT &input, ArrayT &output, const ArrayT &sigma_arr) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif

    static Array<double> temp;
    ASSERT2(sigma_arr.sizes == input.sizes, "sigma_arr sizes != input sizes");
    output.resize(input.sizes);
    cumsum_xy<ArrayT, real, atom, channels, Array<double>, order>(input, temp);
//    for (int i = 1; i < order; i++) {
//        cumsum_xy<Array<double>, real, atom, channels, Array<double> >(temp, temp);
//    }
//    double T1 = wall_time();
//    save_color_image(temp, "gaussian_blur_sat_cumsum.pfm");

    int w = output.width(), h = output.height();
    #pragma omp parallel for
    for (int y = 0; y < h; y++) {
        if (order == 1) {
            for (int x = 0; x < w; x++) {
                GET_SPACING_SCALE();
                int y1 = y-spacing;
                int y2 = y+spacing;
                int x1 = x-spacing;
                int x2 = x+spacing;

                /* TODO: Should clamp x1, y1, x2, y2 first and then use get_nearest(). */
                for (int c = 0; c < channels; c++) {
                    output.get_nearest(y, x, c) = (temp.get_clamp(y1,x1,c) - temp.get_clamp(y1,x2,c) - temp.get_clamp(y2,x1,c) + temp.get_clamp(y2,x2,c))*scale;
                }
            }
        } else if (order == 2) {
            for (int x = 0; x < w; x++) {
                GET_SPACING_SCALE();
                
                int y2 = y-1;
                int y1 = y2-spacing;
                int y3 = y2+spacing;

                int x2 = x-1;
                int x1 = x2-spacing;
                int x3 = x2+spacing;

                if (y1 >= 0 && y3 < h && x1 >= 0 && x3 < w) {
                    OUTPUT_ORDER2(0, 0, 0, 0, 0);
                } else {
                    if (y1 < 0) { y1 = 0; if (y2 < 0) { y2 = 0; } }
                    int dy = 0;
                    if (y3 >= h) {
                        dy = y3-(h-1);
                        y3 = h-1;
                    }

                    if (x1 < 0) { x1 = 0; if (x2 < 0) { x2 = 0; } }
                    int dx = 0;
                    if (x3 >= w) {
                        dx = x3-(w-1);
                        x3 = w-1;
                    }
                
                    if (!dy) {
                        OUTPUT_ORDER2((temp.get_nearest(y1, w-1, c)-temp.get_nearest(y1, w-2, c))*dx,
                                      (temp.get_nearest(y2, w-1, c)-temp.get_nearest(y2, w-2, c))*dx,
                                      (temp.get_nearest(y3, w-1, c)-temp.get_nearest(y3, w-2, c))*dx, 0, 0);
                    } else if (!dx) {
                        OUTPUT_ORDER2(0,
                                      0,
                                      (temp.get_nearest(h-1, x3, c)-temp.get_nearest(h-2, x3, c))*dy,
                                      (temp.get_nearest(h-1, x1, c)-temp.get_nearest(h-2, x1, c))*dy,
                                      (temp.get_nearest(h-1, x2, c)-temp.get_nearest(h-2, x2, c))*dy);
                    } else {
                        OUTPUT_ORDER2((temp.get_nearest(y1, w-1, c)-temp.get_nearest(y1, w-2, c))*dx,
                                      (temp.get_nearest(y2, w-1, c)-temp.get_nearest(y2, w-2, c))*dx,
                                      (temp.get_nearest(y3, w-1, c)-temp.get_nearest(y3, w-2, c))*dx+
                                      ((temp.get_nearest(y3,   w-1, c) + (temp.get_nearest(y3,   w-1, c)-temp.get_nearest(y3,   w-2, c))*dx)-
                                       (temp.get_nearest(y3-1, w-1, c) + (temp.get_nearest(y3-1, w-1, c)-temp.get_nearest(y3-1, w-2, c))*dx)
                                      )*dy,
                                      (temp.get_nearest(h-1, x1, c)-temp.get_nearest(h-2, x1, c))*dy,
                                      (temp.get_nearest(h-1, x2, c)-temp.get_nearest(h-2, x2, c))*dy);
                        
                        // Simplification of the third argument to OUTPUT_ORDER2() (not any faster than the original):
                        
                        // a   b         (y3-1,w-2)              (y3-1, w-1)
                        // c   d     =   (y3,  w-2)              (y3,   w-1)
                        
                        // (d-c)*dx + ((d + (d-c)*dx)-(b + (b-a)*dx))*dy =
                        // a*(dx*dy) + b*(-(1+dx)*dy) + c*(-dx*(dy+1)) + d*((1+dx)*dy+dx)
                        
    //                                  (temp.get_nearest(y3, w-1, c)-temp.get_nearest(y3, w-2, c))*dx+
    //                                  ((temp.get_nearest(y3,   w-1, c) + (temp.get_nearest(y3,   w-1, c)-temp.get_nearest(y3,   w-2, c))*dx)-
    //                                   (temp.get_nearest(y3-1, w-1, c) + (temp.get_nearest(y3-1, w-1, c)-temp.get_nearest(y3-1, w-2, c))*dx)
    //                                  )*dy,
    //                         =

    //                                  temp.get_nearest(y3-1, w-2, c)*(dx*dy)-
    //                                  temp.get_nearest(y3-1, w-1, c)*((1+dx)*dy)-
    //                                  temp.get_nearest(y3,   w-2, c)*((1+dy)*dx)+
    //                                  temp.get_nearest(y3,   w-1, c)*((1+dx)*dy+dx),
                    }
                }
            }
        } else {
            ASSERT2(false, "order not implemented in gaussian_blur_sat");
        }
    }
#if DEBUG_TIME
    double T1 = wall_time();
    printf("gaussian_blur_sat: %f\n", T1-T0);
#endif
//    save_color_image(output, "gaussian_blur_sat.pfm");
}

template<class ArrayT, class real, class atom, int channels>
void box_blur(const ArrayT &input, ArrayT &output, int n) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif
    ArrayT temp;
    temp.resize(input.sizes);
    output.resize(input.sizes);

    if (n % 2 == 0) { n++; }
    n = MIN(MIN(n, (input.width()-1)*2), (input.height()-1)*2);

    int n2 = n/2;

#define BOX_BLUR_DELTA 4

    int x1 = n2+1;
    int x2 = input.width()-n2;

#define BOX_BLUR_H_NOCLAMP() \
            for (int c = 0; c < channels; c++) { \
                temp.get_nearest(y, x, c) = temp.get_nearest(y, x-1, c) - input.get_clamp(y, x-n2-1, c) + input.get_clamp(y, x+n2, c); \
            }

    /* Blur x into temp */
    #pragma omp parallel for
    for (int y = 0; y < input.height(); y++) {
        for (int c = 0; c < channels; c++) {
            temp.get_nearest(y, 0, c) = input.get_nearest(y, 0, c)*(n2+1);
        }
        for (int i = 1; i <= n2; i++) {
            for (int c = 0; c < channels; c++) {
                temp.get_nearest(y, 0, c) += input.get_nearest(y, i, c);
            }
        }
//        for (int x = 1; x < input.width(); x++) {
//            BOX_BLUR_H_NOCLAMP();
//        }

        for (int x = 1; x < x1; x++) {
            BOX_BLUR_H_NOCLAMP();
        }
        real *temp_y_0 = temp.data + y*input.stride[0];
        real *temp_y_m_1 = temp.data + y*input.stride[0] + (-1)*input.stride[1];
        real *input_y_m_n2 = input.data + y*input.stride[0] + (-n2-1)*input.stride[1];
        real *input_y_n2 = input.data + y*input.stride[0] + (n2)*input.stride[1];
        for (int x = x1; x < x2; x++) {
            for (int c = 0; c < channels; c++) {
                int offset = x*input.stride[1]+c;
                temp_y_0[offset] = temp_y_m_1[offset] - input_y_m_n2[offset] + input_y_n2[offset];
//                temp.get_nearest(y, x, c) = temp.get_nearest(y, x-1, c) - input.get_nearest(y, x-n2-1, c) + input.get_nearest(y, x+n2, c);
            }
        }
        for (int x = x2; x < input.width(); x++) {
            BOX_BLUR_H_NOCLAMP();
        }
    }
#if DEBUG_TIME
    double Tmid = wall_time();
#endif
#if DEBUG_SAVE_IMAGES
    static int count = 0;
    char buf[256];
    sprintf(buf, "box_blur_mid%d.pfm", count);
    save_color_image(output, buf);
#endif

    /* Blur temp into output */
    int y1 = n2+1;
    int y2 = input.height()-n2;
    
#define BOX_BLUR_V_NOCLAMP(dx_max) \
                for (int dx = 0; dx < dx_max; dx++) { \
                    for (int c = 0; c < channels; c++) { \
                        output.get_nearest(y, x+dx, c) = output.get_nearest(y-1, x+dx, c) + (temp.get_clamp(y+n2, x+dx, c) - temp.get_clamp(y-n2-1, x+dx, c))*scale; \
                    } \
                }

#define BOX_BLUR_V_FULL(dx_max) \
            for (int dx = 0; dx < dx_max; dx++) { \
                for (int c = 0; c < channels; c++) { \
                    output.get_nearest(0, x+dx, c) = temp.get_nearest(0, x+dx, c)*(n2+1)*scale; \
                } \
            } \
            for (int i = 1; i <= n2; i++) { \
                for (int dx = 0; dx < dx_max; dx++) { \
                    for (int c = 0; c < channels; c++) { \
                        output.get_nearest(0, x+dx, c) += temp.get_nearest(i, x+dx, c)*scale; \
                    } \
                } \
            } \
            \
            for (int y = 1; y < y1; y++) { \
                BOX_BLUR_V_NOCLAMP(dx_max); \
            } \
            for (int y = y1; y < y2; y++) { \
                real *output_y = output.data + y*output.stride[0]; \
                real *output_y1 = output.data + (y-1)*output.stride[0]; \
                real *temp_y_n2 = temp.data + (y+n2)*temp.stride[0]; \
                real *temp_y_n1 = temp.data + (y-n2-1)*temp.stride[0]; \
                for (int dx = 0; dx < dx_max; dx++) { \
                    int dx_offset = (x+dx)*output.stride[1]; \
                    for (int c = 0; c < channels; c++) { \
                        int offset = dx_offset + c*output.stride[2]; \
                        output_y[offset] = output_y1[offset] + (temp_y_n2[offset] - temp_y_n1[offset])*scale; \
                    } \
                } \
            } \
            for (int y = y2; y < input.height(); y++) { \
                BOX_BLUR_V_NOCLAMP(dx_max); \
            }

// Inner loop originally:
//                        output.get_nearest(y, x+dx, c) = output.get_nearest(y-1, x+dx, c) + (temp.get_clamp(y+n2, x+dx, c) - temp.get_clamp(y-n2-1, x+dx, c))*scale; \


    double scale = 1.0/(n*n);
    #pragma omp parallel for
    for (int x = 0; x < input.width(); x+=BOX_BLUR_DELTA) {
        if (x + BOX_BLUR_DELTA <= input.width()) {
            BOX_BLUR_V_FULL(BOX_BLUR_DELTA);
        } else {
            int dx_max = input.width() - x;
            BOX_BLUR_V_FULL(dx_max);
        }
    }
#if DEBUG_TIME
    double T1 = wall_time();
    printf("box_blur: %f (%f H, %f V)\n", T1-T0, Tmid-T0, T1-Tmid);
#endif
#if DEBUG_SAVE_IMAGES
    sprintf(buf, "box_blur%d.pfm", count);
    save_color_image(output, buf);
    count++;
#endif
}

#define EUCLIDEAN_DIST_TRY_INBOUNDS(xsrc, ysrc, dx, dy) \
    { \
        int xp_src = output.get_nearest(ysrc, xsrc, X_COORD) + dx; \
        int yp_src = output.get_nearest(ysrc, xsrc, Y_COORD) + dy; \
        int d = xp_src * xp_src + yp_src * yp_src; \
        if (d < dist) { \
            dist = d; \
            xp = xp_src; \
            yp = yp_src; \
        } \
    }

#define EUCLIDEAN_DIST_TRY(xsrc, ysrc, dx, dy) \
    if ((unsigned) (xsrc) < (unsigned) output.width() && (unsigned) (ysrc) < (unsigned) output.height()) { \
        EUCLIDEAN_DIST_TRY_INBOUNDS(xsrc, ysrc, dx, dy); \
    }

#define EUCLIDEAN_DIST_READ() \
            int xp = output.get_nearest(y, x, X_COORD); \
            int yp = output.get_nearest(y, x, Y_COORD); \
            int dist = xp * xp + yp * yp;

#define EUCLIDEAN_DIST_WRITE() \
            output.get_nearest(y, x, X_COORD) = xp; \
            output.get_nearest(y, x, Y_COORD) = yp;

#define EUCLIDEAN_DIST_WRITE_COORD(xsrc, ysrc) \
    if ((unsigned) xsrc < (unsigned) output.width() && (unsigned) ysrc < (unsigned) output.height() && input.get_nearest(ysrc, xsrc, 0)) { \
        xp = xsrc; \
        yp = ysrc; \
        EUCLIDEAN_DIST_WRITE(); \
        continue; \
    }



#define DIST_TRY_INBOUNDS(xsrc, ysrc, dx, dy) \
    { \
        int v_src = output.get_nearest(ysrc, xsrc); \
        int xp_src = INT_TO_X(v_src); \
        int yp_src = INT_TO_Y(v_src); \
        int d = DIST_2D(x, xp_src, y, yp_src); \
        if (d < dist) { \
            dist = d; \
            xp = xp_src; \
            yp = yp_src; \
        } \
    }

#define DIST_TRY(xsrc, ysrc, dx, dy) \
    if ((unsigned) (xsrc) < (unsigned) output.width() && (unsigned) (ysrc) < (unsigned) output.height()) { \
        DIST_TRY_INBOUNDS(xsrc, ysrc, dx, dy); \
    }

#define DIST_READ() \
            int vp = output.get_nearest(y, x); \
            int xp = INT_TO_X(vp); \
            int yp = INT_TO_Y(vp); \
            int dist = DIST_2D(x, xp, y, yp);

#define DIST_WRITE() \
            output.get_nearest(y, x) = XY_TO_INT(xp, yp);

#define DIST_WRITE_COORD(xsrc, ysrc) \
    if ((unsigned) xsrc < (unsigned) output.width() && (unsigned) ysrc < (unsigned) output.height() && input.get_nearest(ysrc, xsrc, 0)) { \
        xp = xsrc; \
        yp = ysrc; \
        DIST_WRITE(); \
        continue; \
    }


#define DIST_2D(ax, bx, ay, by) ((ax-bx)*(ax-bx) + (ay-by)*(ay-by))

/* Euclidean distance transform (returns x, y coordinates). From Jump Flooding paper, 2006. */
template<class ArrayT>
void euclidean_dist_jumpflood(const ArrayT &input, Array<int> &output) {
    vector<int> output_sizes({input.height(), input.width(), 2});
    output.resize(output_sizes);
    int inf = 30*1000;
    
    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        for (int x = 0; x < output.width(); x++) {
            if (!input.get_nearest(y,x,0)) {
                output.get_nearest(y, x, X_COORD) = -1;
                output.get_nearest(y, x, Y_COORD) = -1;
            } else {
                output.get_nearest(y, x, X_COORD) = x;
                output.get_nearest(y, x, Y_COORD) = y;
            }
        }
    }

    int jump = 1;
    while (jump < input.height() || jump < input.width()) {
        jump *= 2;
    }
    
    while (jump >= 1) {
        #pragma omp parallel for
        for (int y = 0; y < output.height(); y++) {
            for (int x = 0; x < output.width(); x++) {
                int xp = output.get_nearest(y, x, X_COORD);
                int yp = output.get_nearest(y, x, Y_COORD);
                int d = xp >= 0 ? DIST_2D(x, xp, y, yp): inf;
                for (int dy = -1; dy <= 1; dy++) {
                    int ysrc = y+dy*jump;
                    if ((unsigned) ysrc >= (unsigned) input.height()) { continue; }
                    for (int dx = -1; dx <= 1; dx++) {
                        int xsrc = x+dx*jump;
                        if ((unsigned) xsrc >= (unsigned) input.width()) { continue; }
                        
                        int xpp = output.get_nearest(ysrc, xsrc, X_COORD);
                        int ypp = output.get_nearest(ysrc, xsrc, Y_COORD);
                        if (xpp >= 0) {
                            int d2 = DIST_2D(x, xpp, y, ypp);
                            if (d2 < d) {
                                xp = xpp;
                                yp = ypp;
                            }
                        }
                    }
                }
                output.get_nearest(y, x, X_COORD) = xp;
                output.get_nearest(y, x, Y_COORD) = yp;
            }
        }
        jump /= 2;
    }
}

/* Euclidean distance transform (returns x, y coordinates). From Euclidean Distance Mapping by Danielsson, 1980. Finds distance to nonzero elements. */
template<class ArrayT>
void euclidean_dist(const ArrayT &input, Array<int> &output) {
//    double T0 = wall_time();
    vector<int> output_sizes({input.height(), input.width(), 2});
    output.resize(output_sizes);
    
    static Array<int> row_x;
    static vector<int> row_size;
    vector<int> row_x_sizes({input.height(), input.width()});
    row_x.resize(row_x_sizes);
    row_size.resize(output.height());
    for (int y = 0; y < output.height(); y++) {
        row_size[y] = 0;
    }
    
    int inf = 30*1000;

    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        for (int x = 0; x < output.width(); x++) {
            if (!input.get_nearest(y,x,0)) {
                output.get_nearest(y, x, X_COORD) = inf;
                output.get_nearest(y, x, Y_COORD) = inf;
                row_x.get_nearest(y, row_size[y]) = x;
                row_size[y]++;
            } else {
                output.get_nearest(y, x, X_COORD) = 0;
                output.get_nearest(y, x, Y_COORD) = 0;
            }
        }
    }

//    double Ta = wall_time();
    for (int y = 0; y < output.height(); y++) {
        for (int ix = 0; ix < row_size[y]; ix++) {
            int x = row_x.get_nearest(y, ix);
            EUCLIDEAN_DIST_READ();

            EUCLIDEAN_DIST_TRY(x-1, y-1, 1, 1);
            EUCLIDEAN_DIST_TRY(x,   y-1, 0, 1);
            EUCLIDEAN_DIST_TRY(x+1, y-1, 1, 1);
            EUCLIDEAN_DIST_TRY(x-1, y,   1, 0);
            
            EUCLIDEAN_DIST_WRITE();
        }
        for (int ix = row_size[y]-1; ix >= 0; ix--) {
            int x = row_x.get_nearest(y, ix);
            EUCLIDEAN_DIST_READ();
            EUCLIDEAN_DIST_TRY(x+1, y,   1, 0);

            EUCLIDEAN_DIST_WRITE();
        }
    }

    for (int y = output.height()-1; y >= 0; y--) {
        for (int ix = 0; ix < row_size[y]; ix++) {
            int x = row_x.get_nearest(y, ix);
            EUCLIDEAN_DIST_READ();
            EUCLIDEAN_DIST_TRY(x-1, y+1, 1, 1);
            EUCLIDEAN_DIST_TRY(x,   y+1, 0, 1);
            EUCLIDEAN_DIST_TRY(x+1, y+1, 1, 1);
            EUCLIDEAN_DIST_TRY(x-1, y,   1, 0);
            
            EUCLIDEAN_DIST_WRITE();
        }

        for (int ix = row_size[y]-1; ix >= 0; ix--) {
            int x = row_x.get_nearest(y, ix);
            EUCLIDEAN_DIST_READ();
            EUCLIDEAN_DIST_TRY(x+1, y, 1, 0);
            
            EUCLIDEAN_DIST_WRITE();
        }
    }
//    double Tb = wall_time();

    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        for (int x = 0; x < output.width(); x++) {
            EUCLIDEAN_DIST_READ();

            EUCLIDEAN_DIST_WRITE_COORD(x+xp, y+yp);
            EUCLIDEAN_DIST_WRITE_COORD(x-xp, y+yp);
            EUCLIDEAN_DIST_WRITE_COORD(x+xp, y-yp);
            EUCLIDEAN_DIST_WRITE_COORD(x-xp, y-yp);
            
/*            printf("%d, %d    %d, %d\n", x, y, xp, yp); */
            ASSERT2(false, "euclidean_dist could not find nonzero sample");
        }
    }

/*
    Array<float> output_x;
    Array<float> output_y;
    output_x.resize(input.sizes);
    output_y.resize(input.sizes);
    output_x.clear();
    output_y.clear();
    for (int y2 = 0; y2 < input.height(); y2++) {
        for (int x2 = 0; x2 < input.width(); x2++) {
            output_x.get_nearest(y2, x2, 0) = output.get_nearest(y2, x2, X_COORD);
            output_y.get_nearest(y2, x2, 0) = output.get_nearest(y2, x2, Y_COORD);
        }
    }
    save_color_image(input, "euclidean_dist_error_input.png");
    save_color_image(output_x, "euclidean_dist_error_output_x.pfm");
    save_color_image(output_y, "euclidean_dist_error_output_y.pfm"); */
//    double T1 = wall_time();
//    printf("euclidean_dist in %f seconds (%f)\n", T1-T0, Tb-Ta);
}

#define COPY_NEAREST_THRESH 1e-44

/* Manhattan distance transform (returns x, y coordinates). From Distance Functions on Digital Pictures, Rosenfeld and Pfaltz 1967. */
template<class ArrayT>
void manhattan_dist(const ArrayT &input, Array<int> &output, double thresh=COPY_NEAREST_THRESH) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif
    vector<int> output_sizes({input.height(), input.width()});
    output.resize(output_sizes);
    
    int inf = 30*1000;

    #pragma omp parallel for
    for (int y = 0; y < output.height(); y++) {
        for (int x = 0; x < output.width(); x++) {
            if (input.get_nearest(y,x,0) < thresh) {
                output.get_nearest(y, x) = XY_TO_INT(inf, inf);
            } else {
                output.get_nearest(y, x) = XY_TO_INT(x, y);
            }
        }
    }

//    double Ta = wall_time();
//    #pragma omp parallel for schedule(static, 8)
    for (int y = 0; y < output.height(); y++) {
        if (y == 0) {
            for (int x = 0; x < output.width(); x++) {
                DIST_READ();

                DIST_TRY(x-1, y,   1, 0);
                
                DIST_WRITE();
            }
        } else {
            {
                int x = 0;
                DIST_READ();

                DIST_TRY(x,   y-1, 0, 1);
                
                DIST_WRITE();
            }

            for (int x = 1; x < output.width(); x++) {
                DIST_READ();

                DIST_TRY_INBOUNDS(x,   y-1, 0, 1);
                DIST_TRY_INBOUNDS(x-1, y,   1, 0);
                
                DIST_WRITE();
            }
        }
    }

//    #pragma omp parallel for schedule(static, 8)
    for (int y = output.height()-1; y >= 0; y--) {
        if (y == output.height()-1) {
            for (int x = output.width()-1; x >= 0; x--) {
                DIST_READ();
                DIST_TRY(x,   y+1, 0, 1);
                DIST_TRY(x+1, y, 1, 0);
                
                DIST_WRITE();
            }
        } else {
            {
                int x = output.width()-1;
                DIST_READ();

                DIST_TRY(x,   y+1, 0, 1);
                
                DIST_WRITE();
            }
            for (int x = output.width()-2; x >= 0; x--) {
                DIST_READ();
                DIST_TRY_INBOUNDS(x,   y+1, 0, 1);
                DIST_TRY_INBOUNDS(x+1, y, 1, 0);
                
                DIST_WRITE();
            }
        }
    }
#if DEBUG_TIME
    double T1 = wall_time();
    printf("manhattan_dist in %f secs\n", T1-T0);
#endif
}

template<class ArrayT>
void copy_nearest(ArrayT &I, ArrayT &sample, int copy_sample=0, double thresh=COPY_NEAREST_THRESH) {
#if DEBUG_TIME
    double T0 = wall_time();
#endif
    ASSERT2(I.dimensions() == 3, "I dimensions not 3");
    ASSERT2(sample.dimensions() == 3, "sample dimensions not 3");

    bool needed = false;
#pragma omp parallel for
    for (int y = 0; y < I.height(); y++) {
        if (!needed) {
            for (int x = 0; x < I.width(); x++) {
                if (sample(y, x, 0) < thresh) {
                    needed = true;
                    break;
                }
            }
        }
    }
    if (!needed) { return; }

    static Array<int> nearest;
    nearest.resize(I.sizes);
    manhattan_dist<ArrayT>(sample, nearest);
    
    if (copy_sample) {
        #pragma omp parallel for
        for (int y = 0; y < I.height(); y++) {
            for (int x = 0; x < I.width(); x++) {
                if (sample(y, x, 0) < thresh) {
                    int v = nearest.get_nearest(y, x);
                    int xsrc = INT_TO_X(v);
                    int ysrc = INT_TO_Y(v);
#if DEBUG
                    if (!in_bounds(xsrc, I.width()) || !in_bounds(ysrc, I.height())) {
                        printf("xsrc,ysrc out of bounds in I: %dx%d\n", xsrc, ysrc); exit(1);
                    }
#endif
                    COPY_FROM(I, xsrc, ysrc);
                    COPY_FROM(sample, xsrc, ysrc);
                }
            }
        }
    } else {
        #pragma omp parallel for
        for (int y = 0; y < I.height(); y++) {
            for (int x = 0; x < I.width(); x++) {
                if (sample(y, x, 0) < thresh) {
                    int v = nearest.get_nearest(y, x);
                    int xsrc = INT_TO_X(v);
                    int ysrc = INT_TO_Y(v);
                    COPY_FROM(I, xsrc, ysrc);
                }
            }
        }
    }
#if DEBUG_TIME
    double T1 = wall_time();
    printf("copy_nearest in %f seconds\n", T1-T0);
#endif
}

template<class T, int T_vec, int channels>
void iir_filter(Array<T, T_vec> &in, const vector<float> &coeffs, float initial=1.0) {
    ASSERT2(coeffs.size() >= 4, "expected at least 4 coefficients to iir_filter()");
    for (int i = 4; i < (int) coeffs.size(); i++) {
        ASSERT2(coeffs[i] == 0.0, "expected at most 4 nonzero coefficients for iir_filter()");
    }
    
    float a0 = coeffs[0];
    float b1 = coeffs[1];
    float b2 = coeffs[2];
    float b3 = coeffs[3];
    int h = in.height(), w = in.width();
    unsigned uh = (unsigned) h, uw = (unsigned) w;
    
#if DEBUG_TIME
    double T0 = wall_time();
#endif

    #pragma omp parallel for schedule(dynamic, 8)
    for (int y = 0; y < h; y += 1) {
        /* Right */
        for (int x = 0; x < 3; x += 1) {
            for (int c = 0; c < channels; c++) {
                T ans = a0*in.get_nearest(y, x, c);
                
                ans += b1*(IIR_IN_BOUNDS(1, 0, 1) ? in.get_nearest(y, x-1, c): initial*in.get_clamp(y, x-1, c));
                ans += b2*(IIR_IN_BOUNDS(1, 0, 2) ? in.get_nearest(y, x-2, c): initial*in.get_clamp(y, x-2, c));
                ans += b3*(IIR_IN_BOUNDS(1, 0, 3) ? in.get_nearest(y, x-3, c): initial*in.get_clamp(y, x-3, c));
                
                in.get_nearest(y, x, c) = ans;
            }
        }

        for (int x = 3; x < w; x += 1) {
            for (int c = 0; c < channels; c++) {
                T ans = a0*in.get_nearest(y, x, c);
                
                ans += b1*(in.get_nearest(y, x-1, c));
                ans += b2*(in.get_nearest(y, x-2, c));
                ans += b3*(in.get_nearest(y, x-3, c));
                
                in.get_nearest(y, x, c) = ans;
            }
        }
        
        /* Left */
        for (int x = w-1; x > w-4; x -= 1) {
            for (int c = 0; c < channels; c++) {
                T ans = a0*in.get_nearest(y, x, c);
                
                ans += b1*(IIR_IN_BOUNDS(-1, 0, 1) ? in.get_nearest(y, x+1, c): initial*in.get_clamp(y, x+1, c));
                ans += b2*(IIR_IN_BOUNDS(-1, 0, 2) ? in.get_nearest(y, x+2, c): initial*in.get_clamp(y, x+2, c));
                ans += b3*(IIR_IN_BOUNDS(-1, 0, 3) ? in.get_nearest(y, x+3, c): initial*in.get_clamp(y, x+3, c));

                in.get_nearest(y, x, c) = ans;
            }
        }

        for (int x = w-4; x >= 0; x -= 1) {
            for (int c = 0; c < channels; c++) {
                T ans = a0*in.get_nearest(y, x, c);
                
                ans += b1*(in.get_nearest(y, x+1, c));
                ans += b2*(in.get_nearest(y, x+2, c));
                ans += b3*(in.get_nearest(y, x+3, c));

                in.get_nearest(y, x, c) = ans;
            }
        }
    }

#if DEBUG_TIME
    double T1 = wall_time();
#endif

    const int tile_width = 16;
    int xtiles = (w+tile_width-1)/tile_width;
    
    #pragma omp parallel for schedule(dynamic, 4)
    for (int xtile = 0; xtile < xtiles; xtile++) {
    
        /* Down */        
        int xmin = xtile*tile_width;
        int xmax = MIN((xtile+1)*tile_width, w);
        for (int y = 0; y < h; y += 1) {
            for (int x = xmin; x < xmax; x += 1) {
                for (int c = 0; c < channels; c++) {
                    T ans = a0*in.get_nearest(y, x, c);
                    
                    ans += b1*(IIR_IN_BOUNDS(0, 1, 1) ? in.get_nearest(y-1, x, c): initial*in.get_clamp(y-1, x, c));
                    ans += b2*(IIR_IN_BOUNDS(0, 1, 2) ? in.get_nearest(y-2, x, c): initial*in.get_clamp(y-2, x, c));
                    ans += b3*(IIR_IN_BOUNDS(0, 1, 3) ? in.get_nearest(y-3, x, c): initial*in.get_clamp(y-3, x, c));
                    
                    in.get_nearest(y, x, c) = ans;
                }
            }
        }
        
        /* Up */        
        for (int y = h-1; y >= 0; y -= 1) {
            for (int x = xmax-1; x >= xmin; x -= 1) {
                for (int c = 0; c < channels; c++) {
                    T ans = a0*in.get_nearest(y, x, c);
                    
                    ans += b1*(IIR_IN_BOUNDS(0, -1, 1) ? in.get_nearest(y+1, x, c): initial*in.get_clamp(y+1, x, c));
                    ans += b2*(IIR_IN_BOUNDS(0, -1, 2) ? in.get_nearest(y+2, x, c): initial*in.get_clamp(y+2, x, c));
                    ans += b3*(IIR_IN_BOUNDS(0, -1, 3) ? in.get_nearest(y+3, x, c): initial*in.get_clamp(y+3, x, c));
                    
                    in.get_nearest(y, x, c) = ans;
                }
            }
        }

    }

#if DEBUG_TIME
    double T2 = wall_time();

    printf("iir_filter right/left: %f secs\n", T1-T0);
    printf("iir_filter down/up: %f secs\n", T2-T1);
    printf("iir_filter: %f secs\n", T2-T0);
#endif
}

template<class T, int T_vec=-1, int channels=3>
void gaussian_blur(Array<T, T_vec> &in, double sigma) {
    iir_filter<T, T_vec, channels>(in, gaussian_kernel(sigma));
}

#if __GNUC__
typedef float array_float4 __attribute__ ((vector_size (16)));
#endif

#if ARRAY_USE_OPENCV
template<class real>
Array<real> imresize(const Array<real> &A, int w, int h, int interpolation=cv::INTER_LANCZOS4, bool prefilter_gaussian_opt=true) {
    cv::Mat input;
    cv::Mat output;
    
    if (prefilter_gaussian_opt && (A.height() > h || A.width() > w)) {
        double sigma_y = double(A.height())/double(h);
        double sigma_x = double(A.width())/double(w);
        if (sigma_y < 1) { sigma_y = 0; }
        if (sigma_x < 1) { sigma_x = 0; }
        sigma_x *= 0.5;
        sigma_y *= 0.5;
#if DEBUG_TIME
        double T0 = wall_time();
#endif
        double mean_sigma = (sigma_x+sigma_y)*0.5;
        if (mean_sigma < 8) {
            input = A.to_cv();
            if (sigma_x > 0 || sigma_y > 0) {
                cv::GaussianBlur(input, input, cv::Size(0, 0), sigma_x, sigma_y, cv::BORDER_REPLICATE);
            }
        } else {
            bool handled_vectorized = false;
            if (A.channels() == 3) {
#if __GNUC__
                Array<array_float4, 2> Acopy;
                Acopy.resize(A.height(), A.width(), 1);
                #pragma omp parallel for
                for (int y = 0; y < A.height(); y++) {
                    for (int x = 0; x < A.width(); x++) {
                        real *A_pixel = A.data + y*A.stride[0] + x*A.stride[1];
                        array_float4 *Acopy_pixel = Acopy.data + y*Acopy.stride[0] + x*Acopy.stride[1];
                        (*Acopy_pixel)[0] = A_pixel[0];
                        (*Acopy_pixel)[1] = A_pixel[1];
                        (*Acopy_pixel)[2] = A_pixel[2];
                        (*Acopy_pixel)[3] = 0.0;
                    }
                }
                gaussian_blur<array_float4, 2, 1>(Acopy, mean_sigma);
                input = Acopy.to_cv_vectorized<real>();
                handled_vectorized = true;
#endif
            }
            
            if (!handled_vectorized) {
                Array<real> Acopy;
                Acopy.assign(A);
                if (Acopy.channels() == 3) {
                    gaussian_blur<real, -1, 3>(Acopy, mean_sigma);
                } else if (Acopy.channels() == 1) {
                    gaussian_blur<real, -1, 1>(Acopy, mean_sigma);
                } else {
                    fprintf(stderr, "gaussian_blur: channels unsupported: %d\n", Acopy.channels());
                }
                input = Acopy.to_cv();
            }
        }
#if DEBUG_TIME
        printf("imresize: prefilter %f\n", wall_time()-T0);
#endif
        //        save_color_image<real>(input, "temp_blurred_opencv.png");
    } else {
        input = A.to_cv();
    }
    
    cv::resize(input, output, cv::Size(w, h), 0.0, 0.0, interpolation);
    return Array<real>(output);
}

template<class real>
Array<real> imdownsample(const Array<real> &A, double factor) {
    return imresize(A, int(A.width()*factor+0.5), int(A.height()*factor+0.5));
}

#endif

#endif
