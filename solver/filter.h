
#ifndef _filter_h
#define _filter_h

#include "array.h"
#include "util.h"
#include "params.h"
#include "timer.h"
#include <memory>
#include <utility>

using std::shared_ptr;

/* ------------------------------------------------------------------------
   Filter: Base class for image filter
   ------------------------------------------------------------------------ */

using std::static_pointer_cast;
using std::dynamic_pointer_cast;
using std::to_string;
using std::pair;

typedef vector<shared_ptr<Array<int> > > boxes_t;

enum MaskType {
    MASK_ZERO,
    MASK_VAR,
    MASK_CONSTANT
};

template<class real>
class Filter { public:
    virtual vector<real *> coeffs(bool remove_mask=true) = 0;        /* Vector of pointer to each parameter. Coeffs can be removed due to masking (see FilterSparse) */
    virtual int ncoeffs(bool remove_mask=true) { return coeffs(remove_mask).size(); }     /* Optional, number of coefficients */
    virtual boxes_t coeffs_boxes() = 0;         /* Vector of 1D or 2D arrays of parameters, containing ints storing the index to each coeff. All coeffs are reported, even those that are masked (see FilterSparse). */
    virtual void apply(const vector<const Array<adouble> *> &in, Array<adouble> &out, vector<Array<adouble> > &temp_images) = 0;
    virtual void apply(const vector<const Array<double> *> &in, Array<double> &out, vector<Array<double> > &temp_images) = 0;
    virtual string str(bool flatten=false) = 0;
    virtual shared_ptr<Filter> copy() = 0;
    virtual void sparsity_changed(const vector<MaskType> &sparsity, int offset, int recurse_level) { }
    
    virtual bool is_dag() { return false; }
    virtual bool is_downsample() { return false; }
    virtual bool is_upsample() { return false; }
    virtual bool is_fir() { return false; }
    virtual bool is_iir() { return false; }
    virtual bool is_add() { return false; }
    
    virtual ~Filter() { }
    
    virtual double feature_passes() { return 1; }
    
    virtual double feature_taps(int args) {
        vector<real *> c = coeffs();
        int ans = 0;
        for (int i = 0; i < (int) c.size(); i++) {
            if (to_double(*c[i]) != 0) {
                ans++;
            }
        }
        return ans;
    }

    virtual void feature_surface(const vector<real *> c, const boxes_t &boxes, double &surfacex, double &surfacey) {
        surfacex = 0;
        surfacey = 0;
        for (int i = 0; i < (int) boxes.size(); i++) {
            ASSERT(boxes[i]->dimensions() == 2, "expected 2D box");
            for (int y = 0; y < boxes[i]->height(); y++) {
                for (int x = 0; x < boxes[i]->width(); x++) {
                    int idx = (*boxes[i])(y,x);
                    ASSERT((unsigned) idx < (int) c.size(), "idx out of bounds in feature_surfacex()");
                    if (to_double(*c[idx]) != 0) {
                        if (x > 0) {
                            ASSERT((unsigned) (*boxes[i])(y,x-1) < c.size(), "idx(y,x-1) out of bounds in feature_surface()");
                            if (to_double(*c[(*boxes[i])(y,x-1)]) == 0) { surfacex++; }
                        } else {
                            surfacex++;
                        }
                        if (x+1 < boxes[i]->width()) {
                            ASSERT((unsigned) (*boxes[i])(y,x+1) < c.size(), "idx(y,x+1) out of bounds in feature_surface()");
                            if (to_double(*c[(*boxes[i])(y,x+1)]) == 0) { surfacex++; }
                        } else {
                            surfacex++;
                        }
                        if (y > 0) {
                            ASSERT((unsigned) (*boxes[i])(y-1,x) < c.size(), "idx(y-1,x) out of bounds in feature_surface()");
                            if (to_double(*c[(*boxes[i])(y-1,x)]) == 0) { surfacey++; }
                        } else {
                            surfacey++;
                        }
                        if (y+1 < boxes[i]->height()) {
                            ASSERT((unsigned) (*boxes[i])(y+1,x) < c.size(), "idx(y+1,x) out of bounds in feature_surface()");
                            if (to_double(*c[(*boxes[i])(y+1,x)]) == 0) { surfacey++; }
                        } else {
                            surfacey++;
                        }
                    }
                }
            }
        }
    }

    virtual void feature_stddev(const vector<real *> c, const boxes_t &boxes, double &stddevx, double &stddevy) {
        stddevx = 0;
        stddevy = 0;
        for (int i = 0; i < (int) boxes.size(); i++) {
            ASSERT(boxes[i]->dimensions() == 2, "expected 2D box");
            double meanx_sum = 0.0, meany_sum = 0.0;
            int count = 0;
            for (int y = 0; y < boxes[i]->height(); y++) {
                for (int x = 0; x < boxes[i]->width(); x++) {
                    int idx = (*boxes[i])(y,x);
                    ASSERT((unsigned) idx < (int) c.size(), "idx out of bounds in feature_stddev()");
                    if (to_double(*c[idx]) != 0) {
                        meanx_sum += x;
                        meany_sum += y;
                        count++;
                    }
                }
            }
            if (!count) {
                continue;
            }
            meanx_sum /= count;
            meany_sum /= count;
            
            double currentx = 0.0;
            double currenty = 0.0;
            for (int y = 0; y < boxes[i]->height(); y++) {
                for (int x = 0; x < boxes[i]->width(); x++) {
                    int idx = (*boxes[i])(y,x);
                    ASSERT((unsigned) idx < (int) c.size(), "idx out of bounds in feature_stddev()");
                    if (to_double(*c[idx]) != 0) {
                        double dx = (x-meanx_sum);
                        double dy = (y-meany_sum);
                        currentx += dx*dx;
                        currenty += dy*dy;
                    }
                }
            }
            
            currentx = sqrt(currentx / count);
            currenty = sqrt(currenty / count);

            stddevx += currentx;
            stddevy += currenty;
        }
    }

    virtual double feature_xblocks(const vector<real *> c, const boxes_t &boxes, int blocksize) {
        vector<bool> is_block;
        int ans = 0;
        for (int i = 0; i < (int) boxes.size(); i++) {
            ASSERT(boxes[i]->dimensions() == 2, "expected 2D box");
            is_block.resize((boxes[i]->width()+blocksize-1)/blocksize);
            for (int y = 0; y < boxes[i]->height(); y++) {
                for (int x = 0; x < (int) is_block.size(); x++) {
                    is_block[x] = false;
                }
                
                for (int x = 0; x < boxes[i]->width(); x++) {
                    int idx = (*boxes[i])(y,x);
                    ASSERT((unsigned) idx < (int) c.size(), "idx out of bounds in feature_xblocks()");
                    if (to_double(*c[idx]) != 0) {
                        is_block[x/blocksize] = true;
                    }
                }

                for (int x = 0; x < (int) is_block.size(); x++) {
                    if (is_block[x]) { ans++; }
                }
            }
        }
        return ans;
    }

    virtual vector<double> features(int args=1) {
        double surfacex = 0.0, surfacey = 0.0;

        vector<real *> c = coeffs(false);
        boxes_t boxes = coeffs_boxes();

        feature_surface(c, boxes, surfacex, surfacey);
        if (params.feature == 0) {
            return { feature_passes(), feature_taps(args), surfacex, surfacey };
        } else if (params.feature == 1 || params.feature == 2) {
            const int num_node_types = 5;
            int base_features = 4;
            if (params.feature == 2) { base_features += 5; }
            vector<double> ans(base_features*num_node_types);
            int idx = 0;
            if (is_fir()) { idx = 0; }
            else if (is_iir()) { idx = 1; }
            else if (is_downsample()) { idx = 2; }
            else if (is_upsample()) { idx = 3; }
            else if (is_add()) { idx = 4; }
            else { fprintf(stderr, "node of unknown type: %s\n", str().c_str()); ASSERT2(false, "node of unknown type"); }
            idx *= base_features;
            ans[idx+0] = feature_passes();
            ans[idx+1] = feature_taps(args);
            ans[idx+2] = surfacex;
            ans[idx+3] = surfacey;
            if (params.feature == 2) {
                double stddevx = 0.0, stddevy = 0.0;
                feature_stddev(c, boxes, stddevx, stddevy);
                ans[idx+4] = stddevx;
                ans[idx+5] = stddevy;
                
                ans[idx+6] = feature_xblocks(c, boxes, 2);
                ans[idx+7] = feature_xblocks(c, boxes, 4);
                ans[idx+8] = feature_xblocks(c, boxes, 8);
            }
            return ans;
        } else {
            ASSERT2(false, "feature should be 0 or 1");
        }
    }
    
    virtual double time(int args=1) {
        vector<double> featureL = features(args);
        vector<double> *coeffL = params.feature_coeff_list();
        ASSERT2(featureL.size() == coeffL->size(), "expected featureL and coeffL to be the same size");
        double ans = 0.0;
        for (int i = 0; i < (int) featureL.size(); i++) {
            ans += featureL[i] * (*coeffL)[i];
        }
        return ans;
    }
};

#define PARSE_SINGLE_INPUT(errorMessage) \
ASSERT(in0.size() == 1, errorMessage); \
const Array<T> &in = *in0[0];

#define DECLARE_COPY_SAMECLASS_FUNC(cls, copyfunc) \
shared_ptr<cls > copy_sameclass() { \
    return static_pointer_cast<cls >(copyfunc()); \
}

#define DECLARE_COPY_SAMECLASS(cls) \
DECLARE_COPY_SAMECLASS_FUNC(cls, copy)

#define DECLARE_COPY_APPLY(cls) \
DECLARE_COPY_SAMECLASS_FUNC(cls, copy) \
\
void apply(const vector<const Array<adouble> *> &in0, Array<adouble> &out, vector<Array<adouble> > &temp_images) { \
    apply_type(in0, out, temp_images); \
} \
\
void apply(const vector<const Array<double> *> &in0, Array<double> &out, vector<Array<double> > &temp_images) { \
    apply_type(in0, out, temp_images); \
}

#define DECLARE_NO_COEFFS() \
vector<real *> coeffs(bool remove_mask=true) { \
    return {}; \
} \
int ncoeffs(bool remove_mask=true) { \
    return 0; \
} \
virtual void sparsity_changed(const vector<MaskType> &sparsity, int offset, int recurse_level) { \
} \
boxes_t coeffs_boxes() { \
    return {}; \
}

extern double T_fir;

/* ------------------------------------------------------------------------
   FilterFIR: Finite impulse response filter
   ------------------------------------------------------------------------ */
   
template<class real>
class FilterFIR: public Filter<real> { public:
    Array<real> K;
    vector<pair<int, int> > dxdy_L;
    vector<MaskType> mask;
    
    virtual ~FilterFIR() { }
    FilterFIR(const Array<real> &K_): K(K_) {
    }

    bool is_fir() { return true; }

    string str(bool flatten=false) {
        return string("\"FIRFilter\", \n") + K.str();
    }
    
    vector<real *> coeffs(bool remove_mask=true) {
        vector<real *> ans;
        for (int i = 0; i < K.nelems; i++) {
            ans.push_back(&K.data[i]);
        }
        return ans;
    }
    
    int ncoeffs(bool remove_mask=true) {
        return K.nelems;
    }
    
    boxes_t coeffs_boxes() {
        ASSERT(K.dimensions() == 2, "FilterFIR: expected 2D coeffs");
        auto box = make_shared<Array<int> >(K.sizes);
        for (int i = 0; i < K.nelems; i++) {
            box->data[i] = i;
        }
        return {box};
    }

    void sparsity_changed(const vector<MaskType> &sparsity, int offset, int recurse_level) {
        int nsize = K.nelems;
        //if (offset + nsize > int(sparsity.size())) {
        //    fprintf(stderr, "offset and size out of bounds: %d %d, %d, recurse_level=%d\n", offset, nsize, int(sparsity.size()), recurse_level); ASSERT(false, "");
        //}
        ASSERT(offset + nsize <= int(sparsity.size()), "offset and size out of bounds");
        mask.resize(nsize);
        for (int i = 0; i < nsize; i++) {
            mask[i] = sparsity[offset+i];
        }
    }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        double T0 = wall_time();
        PARSE_SINGLE_INPUT("FilterFIR: expected 1 input");
        //ASSERT(in.dimensions() == 2, "FilterFIR: expected 2D input image");
        ASSERT(in.dimensions() == 2 || in.dimensions() == 3, "FilterFIR: expected 2D or 3D input image");
        ASSERT(K.dimensions() == 2, "FilterFIR: expected 2D coeffs");
//        printf("FilterFIR(%p)::apply_type calls resize\n", this);
        out.resize(in.sizes);
        
        int in_width = in.width();
        int in_height = in.height();
        int hw = K.width()/2;
        int hh = K.height()/2;
        int w = K.width();
        int h = K.height();
                
        int dxdy_count = 0;
        dxdy_L.resize(w*h);
        if (K.nelems != (int) mask.size()) {
            if (mask.size() == 0) {
                mask.resize(K.nelems);
                for (int i = 0; i < K.nelems; i++) { mask[i] = MASK_VAR; }
            } else {
                fprintf(stderr, "Error: K and mask have different number of elements\n");
                fprintf(stderr, "current: %s\n", str().c_str());
                fprintf(stderr, "K elems: %d, mask elems: %d\n", K.nelems, int(mask.size()));
                exit(1);
            }
        }
        ASSERT(K.nelems == (int) mask.size(), "expected K and mask to have same number of elements");
        for (int dy = -hh; dy <= hh; dy++) {
            int Ky = hh-dy;
            if ((unsigned) Ky >= (unsigned) h) { continue; }
            for (int dx = -hw; dx <= hw; dx++) {
                int Kx = (hw-dx);
                if ((unsigned) Kx >= (unsigned) w) { continue; }
                int Kidx = Ky*w+Kx;
                ASSERT((unsigned) Kidx < (unsigned) K.nelems, "Kidx out of bounds");
                if (mask[Kidx] != MASK_ZERO) { dxdy_L[dxdy_count++] = pair<int, int>(dx, dy); }
            }
        }
        
        int numChannels;
        if (in.dimensions() == 2){
            numChannels = 1;
        } else if (in.dimensions() == 3){
            numChannels = out.channels();
        } else {
            ASSERT2(0, "dimensions is not 2 or 3"); exit(1);
        }
        for (int z = 0; z < numChannels; z++) {
            for (int y = 0; y < out.height(); y++) {
                for (int x = 0; x < out.width(); x++) {
                    T ans = 0;
                    for (int i = 0; i < dxdy_count; i++) {
                        int dx = dxdy_L[i].first;
                        int dy = dxdy_L[i].second;
                        int yp = y + dy;
                        int xp = x + dx;
                        if ((unsigned) yp >= (unsigned) in_height ||
                            (unsigned) xp >= (unsigned) in_width) { continue; }
                        if (in.dimensions() == 2){
                            ans += cast_types<real, T>(K(hh-dy, hw-dx)) * in(yp, xp);
                        } else if (in.dimensions() == 3){
                            ans += cast_types<real, T>(K(hh-dy, hw-dx)) * in(yp, xp,z);
                        }
                    }
                    if (in.dimensions() == 2){
                        out(y, x) = ans;
                    } else if (in.dimensions() == 3){
                        out(y, x, z) = ans;
                    }
                }
            }
        }
        T_fir += wall_time() - T0;
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterFIR<real> >(K);
    }

    DECLARE_COPY_APPLY(FilterFIR<real>);
};

template<class real>
class FilterFIRH: public FilterFIR<real> { public:
    FilterFIRH(const Array<real> &K_)
    :FilterFIR<real>(K_) { }

    virtual ~FilterFIRH() { }
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterFIRH<real> >(FilterFIR<real>::K);
    }
    
    DECLARE_COPY_SAMECLASS(FilterFIRH<real>);
};

template<class real>
class FilterFIRV: public FilterFIR<real> { public:
    FilterFIRV(const Array<real> &K_)
    :FilterFIR<real>(K_) { }

    virtual ~FilterFIRV() { }
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterFIRV<real> >(FilterFIR<real>::K);
    }
    
    DECLARE_COPY_SAMECLASS(FilterFIRV<real>);
};

/* ------------------------------------------------------------------------
   FilterAdd: Adds 2+ images
   ------------------------------------------------------------------------ */

template<class real>
class FilterAdd: public Filter<real> { public:
    
    virtual ~FilterAdd() { }
    FilterAdd(){}

    bool is_add() { return true; }

    double feature_taps(int args) {
        return args;
    }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
//        printf("FilterAdd::apply_type calls out.resize\n");
        out.resize(in0[0]->sizes); //ASSERT(in0.size() == 1, filterName+": expected 1 input"); 
//        printf("  FilterAdd::done with out.resize\n");
        out.clear();
        
        for (int i = 0; i < in0.size(); i++) {
			const Array<T> &in = *(in0[i]);
			ASSERT(in.dimensions() == 2, "FilterAdd: expected vector of 2D input images");
			ASSERT(in.sizes == out.sizes, "FilterAdd: expected vector of 2D input images of equal size");
            for (int y = 0; y < out.height(); y++) {
        		for (int x = 0; x < out.width(); x++) {
        			out(y, x) += in(y,x);
        		}
    		}
		}
        //out.normalize2();
	}
	
    DECLARE_NO_COEFFS();
    
    string str(bool flatten=false) {
    	return string("\"FilterAdd\", \n");
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterAdd<real> >();
    }

    DECLARE_COPY_APPLY(FilterAdd<real>);
};

/* ------------------------------------------------------------------------
   FilterTranspose: Image Transpose
   ------------------------------------------------------------------------ */

template<class real>
class FilterTranspose: public Filter<real> { public:
    
    virtual ~FilterTranspose() { }
    FilterTranspose(){}
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        PARSE_SINGLE_INPUT("FilterTranspose: expected 1 input");
        ASSERT(in.dimensions() == 2, "FilterTranspose: expected 2D input image");
        out.resize( {in.sizes[1], in.sizes[0]} );
        
        for (int y = 0; y < out.height(); y++) {
            for (int x = 0; x < out.width(); x++) {
                out(y,x) = in(x,y);
            }
        }
        //out.normalize2();
    }
    
    DECLARE_NO_COEFFS();
    
    string str(bool flatten=false) {
    	return string("\"FilterTranspose\", \n");
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterTranspose<real> >();
    }

    DECLARE_COPY_APPLY(FilterTranspose<real>);
};

/* ------------------------------------------------------------------------
   FilterFlipH: Horizontal Image Flip
   ------------------------------------------------------------------------ */

template<class real>
class FilterFlipH: public Filter<real> { public:
    
    virtual ~FilterFlipH() { }
    FilterFlipH(){}
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        PARSE_SINGLE_INPUT("FilterFlipH: expected 1 input");
        ASSERT(in.dimensions() == 2, "FilterFlipH: expected 2D input image");
        out.resize(in.sizes);
        
        for (int y = 0; y < in.height(); y++) {
            for (int x = 0; x < in.width(); x++) {
                out(y, x) = in(y,in.width()-1-x);
            }
        }
        //out.normalize2();
    }
    
    DECLARE_NO_COEFFS();
    
    string str(bool flatten=false) {
    	return string("\"FilterFlipH\", \n");
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterFlipH<real> >();
    }

    DECLARE_COPY_APPLY(FilterFlipH<real>);
};

/* ------------------------------------------------------------------------
   FilterUpsample: Upsample
   ------------------------------------------------------------------------ */

template<class real>
class FilterUpsample: public Filter<real> { public:
    
    virtual ~FilterUpsample() { }
    FilterUpsample(){}
    
    bool is_upsample() { return true; }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        PARSE_SINGLE_INPUT("FilterUpsample: expected 1 input");
        ASSERT(in.dimensions() == 2, "FilterUpsample: expected 2D input image");
//        printf("FilterUpsample::apply_type calls resize\n");
        if (out.dimensions() != 2 || out.height() != in.height()*2 || out.width() != in.width()*2) {
            out.resize( {in.sizes[0] * 2, in.sizes[1] * 2} );
        }
        
        for (int y = 0; y < in.height(); y++) {
            for (int x = 0; x < in.width(); x++) {
                out(y*2, x*2) = in(y,x);
                out(y*2, x*2+1) = 0.0;
                out(y*2+1, x*2) = 0.0;
                out(y*2+1, x*2+1) = 0.0;
            }
        }
        //out.normalize2();
    }
    
    DECLARE_NO_COEFFS();
    
    string str(bool flatten=false) {
    	return string("\"FilterUpsample\", \n");
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterUpsample<real> >();
    }
	
    DECLARE_COPY_APPLY(FilterUpsample<real>);
};

/* ------------------------------------------------------------------------
   FilterDownsample: Downsample
   ------------------------------------------------------------------------ */

template<class real>
class FilterDownsample: public Filter<real> { public:
    
    virtual ~FilterDownsample() { }
    FilterDownsample(){}

    bool is_downsample() { return true; }

    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        PARSE_SINGLE_INPUT("FilterDownsample: expected 1 input");
        ASSERT2(in.dimensions() == 2, "FilterDownsample: expected 2D input image");
        ASSERT2( (in.sizes[0] % 2) == 0 && (in.sizes[1] % 2) == 0, "FilterDownsample: expected 2D input image with even dimensions");
//        printf("FilterDownsample::apply_type calls resize\n");
        if (out.dimensions() != 2 || out.width() != in.width()/2 || out.height() != in.height()/2) {
            out.resize( {in.sizes[0] / 2, in.sizes[1] / 2} );
        }
        
        for (int y = 0; y < in.height(); y = y + 2 ) {
            for (int x = 0; x < in.width(); x = x + 2 ) {
                out(y/2, x/2) = in(y,x);
            }
        }
        //out.normalize2();
    }
    
    DECLARE_NO_COEFFS();
    
    string str(bool flatten=false) {
    	return string("\"FilterDownsample\", \n");
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterDownsample<real> >();
    }
	
    DECLARE_COPY_APPLY(FilterDownsample<real>);
};

/* ------------------------------------------------------------------------
   FilterIIR: Infinite impulse response filter
   ------------------------------------------------------------------------ */

template<class real>
class FilterIIR: public Filter<real> { public:
    Array<real> K;  // Coeffs for FIR filter
    Array<real> F;  // Coeffs for IIR filter (feedback coeffs)
    int xdir, ydir; // Direction of filter, one of (+1, 0), (-1, 0), (0, +1), (0, -1)
    
    virtual ~FilterIIR() { }

    bool is_iir() { return true; }

    FilterIIR(const Array<real> &K_, const Array<real> &F_, int xdir_=1, int ydir_=0): K(K_), F(F_), xdir(xdir_), ydir(ydir_) {
    }
    
    string str(bool flatten=false) {
        return string("\"IIRFilter\", \n") + K.str() + ", " + F.str() + ", " + to_string(xdir) + ", " + to_string(ydir);
    }
    
    vector<real *> coeffs(bool remove_mask=true) {
        vector<real *> ans;
        for (int i = 0; i < K.nelems; i++) {
            ans.push_back(&K.data[i]);
        } 
        for (int i = 0; i < F.nelems; i++) {
            ans.push_back(&F.data[i]);
        } 
        return ans;
    }

    int ncoeffs(bool remove_mask=true) {
        return K.nelems + F.nelems;
    }
    
    boxes_t coeffs_boxes() {
        ASSERT(K.dimensions() == 1, "FilterIIR: expected 1D feed-forward coeffs");
        ASSERT(F.dimensions() == 1, "FilterIIR: expected 1D feed-back coeffs");
        vector<int> boxK_sizes({ydir ? K.nelems: 1, ydir ? 1: K.nelems});
        vector<int> boxF_sizes({ydir ? F.nelems: 1, ydir ? 1: F.nelems});
        auto boxK = make_shared<Array<int> >(boxK_sizes);
        auto boxF = make_shared<Array<int> >(boxF_sizes);
        for (int i = 0; i < K.nelems; i++) {
            boxK->data[i] = i;
        }
        for (int i = 0; i < F.nelems; i++) {
            boxF->data[i] = i + K.nelems;
        }
        return {boxK, boxF};
    }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        PARSE_SINGLE_INPUT("FilterIIR: expected 1 input");
        ASSERT(in.dimensions() == 2, "FilterIIR: expected 2D input image");
        ASSERT(K.dimensions() == 1, "FilterIIR: expected 1D feed-forward coeffs");
        ASSERT(F.dimensions() == 1, "FilterIIR: expected 1D feed-back coeffs");
//        printf("FilterIIR::apply_type calls resize\n");
        out.resize(in.sizes);
        
        //ASSERT (((ydir == 1 || ydir == -1) && xdir == 0) || ((xdir == 1 || xdir == -1) && ydir == 0), "xdir, ydir must be (1, 0), (-1, 0), (0, 1), (0, -1)");
        int ystart = 0, yend = out.height(), ystep = 1;
        int xstart = 0, xend = out.width(), xstep = 1;
        if (ydir < 0) { ystart = out.height()-1; yend = -1; ystep = -1; }
        if (xdir < 0) { xstart = out.width()-1; xend = -1; xstep = -1; }
        int hw = K.size()/2;
        int w = K.size();
        
        for (int y = ystart; y != yend; y += ystep) {
            for (int x = xstart; x != xend; x += xstep) {
                T ans = 0;
                for (int dx = -hw; dx <= hw; dx++) {
                    //int yp = y-dx*ydir, xp = x-dx*xdir;
                    int yp = y + (ydir ? dx: 0);
                    int xp = x + (xdir ? dx: 0);
                    if ((unsigned) yp < (unsigned) in.height() && (unsigned) xp < (unsigned) in.width() && (unsigned) (hw-dx) < (unsigned) w) {
                        ans += in(yp, xp) * cast_types<real, T>(K(hw-dx));
                    }
                }
                for (int i = 0; i < F.size(); i++) {
                    int yp = y-(i+1)*ydir, xp = x-(i+1)*xdir;
                    if ((unsigned) yp < (unsigned) in.height() && (unsigned) xp < (unsigned) in.width()) {
                        ans += out(y-(i+1)*ydir,x-(i+1)*xdir) * cast_types<real, T>(F(i));
                    }
                }
                out(y, x) = ans;
            }
        }
        //out.normalize2();
    }
    
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterIIR<real> >(K,F,xdir,ydir);
    }

    DECLARE_COPY_APPLY(FilterIIR<real>);
};

template<class real>
class FilterIIRMinusX: public FilterIIR<real> { public:
    FilterIIRMinusX(const Array<real> &K_, const Array<real> &F_)
    :FilterIIR<real>(K_, F_, -1, 0) { }
    
    virtual ~FilterIIRMinusX() { }
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterIIRMinusX<real> >(FilterIIR<real>::K, FilterIIR<real>::F);
    }
    
    DECLARE_COPY_SAMECLASS(FilterIIRMinusX<real>);
};

template<class real>
class FilterIIRPlusY: public FilterIIR<real> { public:
    FilterIIRPlusY(const Array<real> &K_, const Array<real> &F_)
    :FilterIIR<real>(K_, F_, 0, 1) { }

    virtual ~FilterIIRPlusY() { }
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterIIRPlusY<real> >(FilterIIR<real>::K, FilterIIR<real>::F);
    }
    
    DECLARE_COPY_SAMECLASS(FilterIIRPlusY<real>);
};

template<class real>
class FilterIIRMinusY: public FilterIIR<real> { public:
    FilterIIRMinusY(const Array<real> &K_, const Array<real> &F_)
    :FilterIIR<real>(K_, F_, 0, -1) { }

    virtual ~FilterIIRMinusY() { }
    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterIIRMinusY<real> >(FilterIIR<real>::K, FilterIIR<real>::F);
    }
    
    DECLARE_COPY_SAMECLASS(FilterIIRMinusY<real>);
};

/* ------------------------------------------------------------------------
   FilterDAG: DAG of filters
   ------------------------------------------------------------------------ */

#define INPUT_SOURCE_IMAGE (-1)                     /* DAG index for the source image */

#define INPUT_SIZE(in_idx) ((in_idx) >= 0 ? input_size[in_idx]: 1.0)

template<class real>
class NodeDAG { public:
    shared_ptr<Filter<real> > f;
    vector<int> inputs;                             /* Indices of inputs to filter in DAG list */
    
    NodeDAG() :f(NULL) { }
    NodeDAG(shared_ptr<Filter<real> > f_, const vector<int> &inputs_) :f(f_), inputs(inputs_) { }
    
    NodeDAG copy() {
        return NodeDAG(f->copy(), inputs);
    }
};

template<class real>
class FilterDownsamplePrefilter;

template<class real>
class FilterUpsamplePostfilter;

template<class real>
class FilterDAG: public Filter<real> { public:
    vector<NodeDAG<real> > L;                              /* Sequence of filtered images. Last is output. */

    vector<vector<Array<double> > >  temp_images_subL_double;
    vector<vector<Array<adouble> > > temp_images_subL_adouble;

    virtual ~FilterDAG() { }
    
    FilterDAG() { }
    FilterDAG(const vector<NodeDAG<real> > &L_, bool ignore) :L(L_) { }
    FilterDAG(const vector<shared_ptr<Filter<real> > > &L_) {
        L.resize(L_.size());
        for (int i = 0; i < (int) L_.size(); i++) {
            L[i] = NodeDAG<real>(L_[i], {i-1});
        }
    }
    
    bool is_dag() {
        return true;
    }
    
    void check() {
        //ASSERT(L.size() >= 1, "FilterDAG::check(): contains 0 nodes");
        for (int i = 0; i < (int) L.size(); i++) {
            for (int j = 0; j < (int) L[i].inputs.size(); j++) {
                int idx = L[i].inputs[j];
                int max_idx = MIN((int) L.size(), i);
                if (idx != INPUT_SOURCE_IMAGE && (idx < 0 || idx >= max_idx)) {
                    fprintf(stderr, "FilterDAG: check failed: L[%d].inputs[%d]=%d, out of range %d (%d), for DAG:\n", i, j, idx, (int) L.size(), max_idx);
                    fprintf(stderr, "%s\n\n", str().c_str());
                    exit(1);
                }
            }
        }
    }
    
    /* Check that input/output sizes are compatible */
    bool check_sizes(vector<double> &input_size) {
        const double epsilon = 1e-8;
        input_size.resize(L.size());
        ASSERT(L.size(), "expected L size nonzero in check_sizes()");
        for (int i = 0; i < (int) L.size(); i++) {
            ASSERT(L[i].inputs.size(), "Expected nonzero inputs size");
            int in_idx = L[i].inputs[0];
            ASSERT(in_idx < i, "Expected input idx < i");
            double current_size = INPUT_SIZE(in_idx);
            
            if (L[i].f->is_upsample()) {          /* Also disallow downsample immediately followed by upsample */
                for (int j = 0; j < (int) L[i].inputs.size(); j++) {
                    int input = L[i].inputs[j];
                    if (input >= 0 && L[input].f->is_downsample()) {
                        return false;
                    }
                }
            }
            for (int j = 1; j < (int) L[i].inputs.size(); j++) {
                double size_p = INPUT_SIZE(L[i].inputs[j]);
                if (fabs(current_size - size_p) > epsilon) {
                    return false;
                }
            }
            if (L[i].f->is_downsample()) {
                current_size *= 0.25;
            } else if (L[i].f->is_upsample()) {
                current_size *= 4;
            }
            if (current_size > 1+epsilon) {
                return false;
            }
            input_size[i] = current_size;
        }
        return fabs(input_size[input_size.size()-1]-1) < epsilon;
    }
    
    string index_to_name(int i) {
        if (i == INPUT_SOURCE_IMAGE) {
            return "Input";
        } else if (i == L.size()-1) {
            return "OUT";
        } else {
            return "Filter" + to_string(i);
        }
    }

    void sparsity_changed(const vector<MaskType> &sparsity, int offset, int recurse_level) {
        for (int i = 0; i < (int) L.size(); i++) {
            int ncoeffs = L[i].f->ncoeffs(false);
            //printf("FilterDAG::sparsity_changed: sparsity length: %d, offset: %d, i=%d, ncoeffs=%d, recurse_level=%d\n", int(sparsity.size()), offset, i, ncoeffs, recurse_level);
            L[i].f->sparsity_changed(sparsity, offset, recurse_level+1);
            offset += ncoeffs;
        }
    }

    /* Flatten FilterDAG which contains FilterDAG */
    shared_ptr<FilterDAG<real> > flatten() {
        auto ans = make_shared<FilterDAG<real> >();
        vector<int> old_to_new(L.size(), -1);                 /* Index mapping */
        int out_size = 0;
        for (int i = 0; i < (int) L.size(); i++) {
            if (L[i].f->is_dag()) {
                shared_ptr<FilterDAG<real> > current = static_pointer_cast<FilterDAG<real> >(L[i].f);
                ASSERT(L[i].inputs.size() == 1, "expected 1 input to FilterDAG");
                int input = L[i].inputs[0];
                if (input >= 0) {
                    ASSERT(input < i, "expected input < i");
//                    int input0 = input;
                    input = old_to_new[input];
//                    if ((unsigned) input >= ans->L.size()) {
//                        printf("input out of range: i=%d, input0=%d, input=%d, ans->L.size: %d\n", i, input0, input, (int) ans->L.size());
//                    }
                    ASSERT((unsigned) input < ans->L.size(), "expected input in range");
                }
                for (int j = 0; j < (int) current->L.size(); j++) {
                    ans->L.push_back(current->L[j]);
                    NodeDAG<real> *node = &ans->L[ans->L.size()-1];
                    for (int k = 0; k < (int) node->inputs.size(); k++) {
                        if (node->inputs[k] < 0) {
                            node->inputs[k] = input;
                        } else {
                            node->inputs[k] += out_size;
                        }
                    }
                }
            } else {
                ans->L.push_back(L[i]);
                NodeDAG<real> *node = &ans->L[ans->L.size()-1];
                for (int k = 0; k < (int) node->inputs.size(); k++) {
                    if (node->inputs[k] >= 0) {
                        ASSERT(node->inputs[k] < i, "expected input < i");
                        node->inputs[k] = old_to_new[node->inputs[k]];
                        ASSERT((unsigned) node->inputs[k] < ans->L.size(), "expected input in range");
                    }
                }
            }
            old_to_new[i] = ans->L.size()-1;
            out_size = ans->L.size();
        }
#if DEBUG
        check();
#endif
        return ans;
    }
    
    bool needs_flatten() {
        for (int i = 0; i < (int) L.size(); i++) {
            if (L[i].f->is_dag()) {
                return true;
            }
        }
        return false;
    }

    virtual vector<double> features(int args=1) {
        vector<double> ans;
//        printf("features() of DAG:\n%s\n\n", str(true).c_str());
        
        if (needs_flatten()) { return flatten()->features(args); }
        vector<double> input_size(L.size());
        if (!check_sizes(input_size)) { ASSERT(false, "check failed in FilterDAG::time()"); }
        for (int i = 0; i < (int) L.size(); i++) {
            vector<double> sub = L[i].f->features((int) L[i].inputs.size());
//            printf("sub[%d] = %s, input_size=%f\n", i, vector_to_str_real(sub).c_str(), input_size[i]);
            for (int j = 0; j < (int) sub.size(); j++) {
                sub[j] *= input_size[i];
            }
            if (!ans.size()) { ans = sub; }
            else {
                ASSERT(sub.size() == ans.size(), "feature sizes should all be equivalent");
                for (int j = 0; j < (int) ans.size(); j++) {
                    ans[j] += sub[j];
                }
            }
//            ans +=  * input_size[i];
        }
//        printf("ans = %s\n", vector_to_str_real(ans).c_str());
        return ans;
        
    }

    virtual double time(int args=1) {
        if (needs_flatten()) { return flatten()->time(); }
        return Filter<real>::time(args);
        /*
        double ans = 0.0;
        vector<double> input_size(L.size());
        if (!check_sizes(input_size)) { ASSERT(false, "check failed in FilterDAG::time()"); }
        for (int i = 0; i < (int) L.size(); i++) {
            ans += L[i].f->time((int) L[i].inputs.size()) * input_size[i];
        }
        return ans;
        */
    }

    virtual string str(bool do_flatten=false) {
        if (do_flatten && needs_flatten()) { return flatten()->str(do_flatten); }
        
        string ans("{\n");
        ans += "    \"Input\": [\"ImageParam\", 0],\n";
        for (int i = 0; i < (int) L.size(); i++) {
            string sub = L[i].f->str(do_flatten);
//            FILE *f0 = fopen("f0.txt", "at");
            string search = "\", ";
//            fprintf(f0, "Original: %s\n", sub.c_str());
            bool is_last = (sub.size() >= 3 && sub.substr(sub.size()-search.size()) == search) ||
                           (sub.size() >= 4 && sub.substr(sub.size()-search.size()-1) == search+"\n");
//            fprintf(f0, "is_last: %d\n", int(is_last));
            string replace = search;
            vector<int> &inputs = L[i].inputs;
            for (int j = 0; j < (int) inputs.size(); j++) {
                replace += string("\"") + index_to_name(inputs[j]);
                if (!(is_last && j == int(inputs.size())-1)) { replace += "\", "; }
                else { replace += "\""; }
            }
//            fprintf(f0, "replace: %s\n", replace.c_str());
            sub = string("[") + str_replace(sub, search, replace) + string("]");
            sub = str_replace(sub, "\n", "\n    ");
//            fprintf(f0, "sub: %s\n\n", sub.c_str());
//            fclose(f0);
            ans += string("    \"") + index_to_name(i) + "\": " + sub;
            
            if (i < (int) L.size()-1) {
                ans += ",\n\n";
            }
        }
        ans += "\n}";
        return ans;
    }
    
    virtual vector<real *> coeffs(bool remove_mask=true) {
        vector<real *> ans;
        for (int i = 0; i < (int) L.size(); i++) {
            vector<real *> sub = L[i].f->coeffs(remove_mask);
            ans.insert(ans.end(), sub.begin(), sub.end());
        }
        return ans;
    }

    virtual int ncoeffs(bool remove_mask=true) {
        int ans = 0;
        for (int i = 0; i < (int) L.size(); i++) {
            ans += L[i].f->ncoeffs(remove_mask);
        }
        return ans;
    }
    
    boxes_t subfilter_coeffs_boxes(int i0, bool accum_all=false) {
        boxes_t ans;
        int n = 0;
        for (int i = 0; i <= (int) i0; i++) {
            boxes_t sub = L[i].f->coeffs_boxes();
            for (int j = 0; j < (int) sub.size(); j++) {
                shared_ptr<Array<int> > box = sub[j];
                for (int k = 0; k < box->nelems; k++) {
                    box->data[k] = n++;
                }
                if (i == i0 || accum_all) {
                    ans.push_back(box);
                }
            }
        }
        return ans;
    }
    
    virtual boxes_t coeffs_boxes() {
        return subfilter_coeffs_boxes(L.size()-1, true);
    }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images, vector<vector<Array<T> > > &temp_images_subL) {
//        printf("FilterDAG::apply_type\n");
        PARSE_SINGLE_INPUT("FilterDAG: expected 1 input");
        if (temp_images.size() != L.size()) {
//            printf("FilterDAG::apply_type calls vector::resize on temp_images (%d to %d)\n", (int) temp_images.size(), (int) L.size());
            temp_images.resize(L.size());
        }
        if (temp_images_subL.size() != L.size()) {
//            printf("FilterDAG::apply_type calls vector::resize on temp_images_subL (%d to %d)\n", (int) temp_images_subL.size() , (int) L.size());
            temp_images_subL.resize(L.size());
        }
        check();
        
//        printf("FilterDAG::apply: before loop\n");
        for (int i = 0; i < (int) L.size(); i++) {
            vector<const Array<T> *> inputs;
            for (int j = 0; j < (int) L[i].inputs.size(); j++) {
                int idx = L[i].inputs[j];
                if (idx == INPUT_SOURCE_IMAGE) {
                    inputs.push_back(&in);
                } else {
                    inputs.push_back(&temp_images[idx]);
                }
            }
            L[i].f->apply(inputs, i == (L.size()-1) ? out: temp_images[i], temp_images_subL[i]);
        }
//        printf("FilterDAG::apply: done\n");
    }
    
    virtual shared_ptr<Filter<real> > copy() {
        vector<NodeDAG<real> > Lcopy;
        Lcopy.reserve(L.size());
        for (int i = 0; i < (int) L.size(); i++) {
            Lcopy.push_back(L[i].copy());
        }
        return make_shared<FilterDAG<real> >(Lcopy, true);
    }

    void apply(const vector<const Array<adouble> *> &in0, Array<adouble> &out, vector<Array<adouble> > &temp_images) {
        apply_type(in0, out, temp_images, temp_images_subL_adouble);
    }

    void apply(const vector<const Array<double> *> &in0, Array<double> &out, vector<Array<double> > &temp_images) {
        apply_type(in0, out, temp_images, temp_images_subL_double);
    }


    DECLARE_COPY_SAMECLASS(FilterDAG<real>);
};

/* ------------------------------------------------------------------------
   Gaussian prefiltering helper functions
   ------------------------------------------------------------------------ */

#define GET_GAUSSIAN_PREFILTER(sz) \
    double sigma = params.prefilter_sigma; \
    Array<real> K; \
    K.resize(sz); \
    real Ksum = 0.0; \
    ASSERT2(params.prefilter_size % 2 == 1, "expected prefilter_size to be odd"); \
    int hw = params.prefilter_size / 2; \
    for (int dx = -hw; dx <= hw; dx++) { \
        int i = dx + hw; \
        K.data[i] = exp(-dx*dx/(2*sigma*sigma)); \
        Ksum += K.data[i]; \
    } \
    for (int i = 0; i < K.nelems; i++) { \
        K.data[i] /= Ksum; \
    }


template<class real>
shared_ptr<Filter<real> > prefilter_gaussian_x() {
    vector<int> sz({1, params.prefilter_size});
    GET_GAUSSIAN_PREFILTER(sz);
    return static_pointer_cast<Filter<real> >(make_shared<FilterFIR<real> >(K));
}

template<class real>
shared_ptr<Filter<real> > prefilter_gaussian_y() {
    vector<int> sz({params.prefilter_size, 1});
    GET_GAUSSIAN_PREFILTER(sz);
    return static_pointer_cast<Filter<real> >(make_shared<FilterFIR<real> >(K));
}

/* ------------------------------------------------------------------------
   FilterDownsamplePrefilter: Downsample after prefilter with 5x5 Gaussian
   ------------------------------------------------------------------------ */

template<class real>
class FilterDownsamplePrefilter: public FilterDAG<real> { public:
    FilterDownsamplePrefilter()
        :FilterDAG<real>({prefilter_gaussian_x<real>(), prefilter_gaussian_y<real>(), shared_ptr<Filter<real> >(new FilterDownsample<real>())})
    { }

    bool is_downsample() { return true; }

    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterDownsamplePrefilter<real> >();
    }

    string str(bool flatten=false) {
    	return string("\"FilterDownsamplePrefilter\", \n");
    }

    DECLARE_NO_COEFFS();
    DECLARE_COPY_SAMECLASS(FilterDownsamplePrefilter<real>);
};

/* ------------------------------------------------------------------------
   FilterUpsamplePostfilter: Upsample then postfilter with 5x5 Gaussian
   ------------------------------------------------------------------------ */

template<class real>
class FilterUpsamplePostfilter: public FilterDAG<real> { public:
    FilterUpsamplePostfilter()
        :FilterDAG<real>({make_shared<FilterUpsample<real> >(), prefilter_gaussian_x<real>(), prefilter_gaussian_y<real>()})
    { }

    bool is_upsample() { return true; }

    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterUpsamplePostfilter<real> >();
    }

    string str(bool flatten=false) {
    	return string("\"FilterUpsamplePostfilter\", \n");
    }

    DECLARE_NO_COEFFS();
    DECLARE_COPY_SAMECLASS(FilterUpsamplePostfilter<real>);
};

int downsample_dim(int w);
Array<int> downsample_box(const Array<int> &box);
void test_multiscale();
vector<int> flatten_boxes(const boxes_t &L);

/* ------------------------------------------------------------------------
   FilterMultires: Wraps any Filter, can solve coefficients in multires
   ------------------------------------------------------------------------ */

#define XY_TO_INT(x, y) ((x)|((y)<<16))
#define INT_TO_X(v) ((v)&(65535))
#define INT_TO_Y(v) ((v)>>16)

boxes_t copy_boxes(const boxes_t &b);

template<class real>
class FilterMultires: public Filter<real> {
public:
    shared_ptr<Filter<real> > filter;
    int level;                          /* Number of times downsampled */
    int nlevels0;

    boxes_t boxes;                      /* Original coeff indices (as a list of boxes, i.e. shared_ptr<Array<int> >) */
    boxes_t d_boxes;                    /* Original coeff indices downsampled to small box size (cluster representatives) */
    boxes_t du_boxes;                   /* Original coeff indices downsampled to small size, then upsampled back up to larger box size */
    vector<int> indices;                /* Corresponding indices */
    vector<int> d_indices;
    vector<int> du_indices;
    
    virtual ~FilterMultires() { }
    
    double feature_passes() {
        return filter->feature_passes();
    }
    
    double feature_taps() {
        copy_coeffs();
        return filter->feature_taps();
    }
    
    void feature_surface(const vector<real *> c, const boxes_t &boxes, double &surfacex, double &surfacey) {
        copy_coeffs();
        filter->feature_surface(c, boxes, surfacex, surfacey);
    }

    /* Copy cluster representative coeffs to all other coeffs */
    void copy_coeffs() {
        vector<real *> f_coeffs(filter->coeffs());
        ASSERT(indices.size() == f_coeffs.size(), "indices size != f_coeffs size");
        if (du_indices.size() != f_coeffs.size()) {
            printf("du_indices size: %d, f_coeffs size: %d\n", (int) du_indices.size(), (int) f_coeffs.size());
        }
        ASSERT(du_indices.size() == f_coeffs.size(), "du_indices size != f_coeffs size");
        for (int i = 0; i < (int) indices.size(); i++) {
            ASSERT(du_indices[i] >= 0 && du_indices[i] < (int) f_coeffs.size(), "du_indices index out of bounds");
            *f_coeffs[i] = *f_coeffs[du_indices[i]];
        }
    }

    vector<real *> coeffs(bool remove_mask=true) {
        vector<real *> ans;
        vector<real *> f_coeffs(filter->coeffs());
        ASSERT(d_indices.size() <= du_indices.size(), "d_indices size > du_indices size");
        for (int i = 0; i < (int) d_indices.size(); i++) {
            ASSERT(d_indices[i] >= 0 && d_indices[i] < int(f_coeffs.size()), "d_indices out of bounds");
            ans.push_back(f_coeffs[d_indices[i]]);
        }
        return ans;
    }

    boxes_t coeffs_boxes() {
        //return copy_boxes(d_boxes);
        boxes_t ans;
        int count = 0;
        for (int i = 0; i < (int) d_boxes.size(); i++) {
            shared_ptr<Array<int> > box(d_boxes[i]);
            shared_ptr<Array<int> > box_ans(make_shared<Array<int> >(box->sizes));
            ASSERT(box->dimensions() == 2, "expected 2D box");
            int box_w = box->width(), box_h = box->height();
            for (int y = 0; y < box_h; y++) {
                for (int x = 0; x < box_w; x++) {
                    (*box_ans)(y, x) = count++;
                }
            }
            ans.push_back(box_ans);
        }
        return ans;
    }
    
    void calc_indices() {
        d_indices = flatten_boxes(d_boxes);
        du_indices = flatten_boxes(du_boxes);
    }
    
    void set_level(int level0, bool init=false) {        /* Set level, 0 is finest scale */
        level = level0;
//        if (!(level >= 0 && level < nlevels0)) {
//            printf("set_level: nlevels=%d, level=%d\n", nlevels0, level);
//        }
//        assert(level >= 0 && level < nlevels0);
        ASSERT(level >= 0 && level < nlevels0, "level out of range");
        if (!init) {
            copy_coeffs();
        }
        d_boxes.clear();
        du_boxes.clear();
        for (int i = 0; i < (int) boxes.size(); i++) {
            Array<int> box = *boxes[i];
            ASSERT(box.dimensions() == 2, "expected 2D box");
            for (int j = 0; j < level; j++) {
                box = downsample_box(box);
            }
            d_boxes.push_back(make_shared<Array<int> >(box));
            du_boxes.push_back(make_shared<Array<int> >(level ? resample_nearest(*d_boxes[i], boxes[i]->sizes, true): *boxes[i]));
        }
        calc_indices();
    }
    
    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
        copy_coeffs();
        filter->apply(in0, out, temp_images);
    }
    
    FilterMultires(const shared_ptr<Filter<real> > &filter_, int level_=0) :filter(filter_) {
//        printf("FilterMultires: wrapping:\n");
//        printf("%s\n", typeid(*filter_).name());
//        printf("%s\n", filter_->str().c_str());
        boxes = filter->coeffs_boxes();

        int max_box = 0;
        for (int i = 0; i < (int) boxes.size(); i++) {
            shared_ptr<Array<int> > box = boxes[i];
            ASSERT(box->dimensions() == 2, "box dimensions not 2");
            max_box = MAX(max_box, box->width());
            max_box = MAX(max_box, box->height());
        }
        
        nlevels0 = 0;
        while (max_box > 1) {
            max_box = downsample_dim(max_box);
            nlevels0++;
        }
        nlevels0 = MAX(nlevels0, 1);
        
        indices = flatten_boxes(boxes);
        set_level(level_, true);
    }

    shared_ptr<Filter<real> > copy() {
        return make_shared<FilterMultires<real> >(filter->copy(), level);
    }
    
    int nlevels() {
        return nlevels0;
    }
    
    string str(bool flatten=false) {
        copy_coeffs();
        return filter->str(flatten);
    }
    
    DECLARE_COPY_APPLY(FilterMultires<real>);
};

/* ------------------------------------------------------------------------
   FilterSparse: Wraps any Filter, masks out coeffs with sparsity mask
   ------------------------------------------------------------------------ */

MaskType int_to_masktype(int x);

template<class real>
class FilterSparse: public Filter<real> {
public:
    vector<MaskType> mask;
    shared_ptr<Filter<real> > filter;
    vector<double> save_coeffs;
    
    virtual ~FilterSparse() { }

    virtual double time(int args=1) {
        shared_ptr<FilterSparse<real> > c = copy_sameclass();
        c->zero_masked_coeffs();
        double ans = c->filter->time();
        return ans;
    }

    virtual vector<double> features(int args=1) {
        shared_ptr<FilterSparse<real> > c = copy_sameclass();
        c->zero_masked_coeffs();
        auto ans = c->filter->features();
        return ans;
    }
   
    void init_mask() {
        mask.resize(filter->coeffs().size(), MASK_VAR);
    }
    
    void mask_from_coeffs() {
        vector<real *> coeffs = filter->coeffs();
        mask.resize(coeffs.size(), MASK_VAR);
        //ASSERT(mask.size() == coeffs.size(), "mask size != filter coeffs size");
        for (int i = 0; i < (int) mask.size(); i++) {
            mask[i] = (*coeffs[i] == 0.0 ? MASK_ZERO: MASK_VAR);
        }
    }
    
    FilterSparse(const shared_ptr<Filter<real> > &filter_, vector<MaskType> *mask_=NULL) :filter(filter_) {
        if (!mask_) {
            init_mask();
        } else {
            mask = *mask_;
        }
    }
    
    virtual vector<real *> coeffs(bool remove_mask=true) {
        vector<real *> ans;
        vector<real *> coeffs = filter->coeffs();
        if (!remove_mask) {
            zero_masked_coeffs();
        }
        ASSERT(mask.size() == coeffs.size(), "mask size != filter coeffs size");
        for (int i = 0; i < (int) mask.size(); i++) {
            if (!remove_mask) {
                ans.push_back(coeffs[i]);
            } else {
                if (mask[i] == MASK_VAR) {
                    ans.push_back(coeffs[i]);
                }
            }
        }
        return ans;
    }
    
    virtual boxes_t coeffs_boxes() {
        return filter->coeffs_boxes();
    }

    template<class T>
    void apply_type(const vector<const Array<T> *> &in0, Array<T> &out, vector<Array<T> > &temp_images) {
//        printf("FilterSparse::apply_type\n");
        vector<real *> coeffs = filter->coeffs();
        save_coeffs.resize(mask.size());
//        printf("FilterSparse::apply: after reserve\n");
        ASSERT(mask.size() == coeffs.size(), "mask size != filter coeffs size");
        for (int i = 0; i < (int) mask.size(); i++) {
            if (mask[i] == MASK_ZERO) {
                save_coeffs[i] = to_double(*coeffs[i]);
                *coeffs[i] = 0;
            }
        }
        //printf("FilterSparse::apply, str=\n%s\n", str(false).c_str());
        //printf("FilterSparse::apply sparsity_changed: mask size: %d\n", int(mask.size()));
        filter->sparsity_changed(mask, 0, 0);
//        printf("FilterSparse::apply: before apply of subfilter\n");
        filter->apply(in0, out, temp_images);
//        printf("FilterSparse::apply: after apply of subfilter\n");
        for (int i = 0; i < (int) mask.size(); i++) {
            if (mask[i] == MASK_ZERO) {
                *coeffs[i] = save_coeffs[i];
            }
        }
//        printf("FilterSparse::apply: done\n");
    }
    
    void zero_masked_coeffs() {
        vector<real *> coeffs = filter->coeffs();
        ASSERT(mask.size() == coeffs.size(), "mask size != filter coeffs size");
        for (int i = 0; i < (int) mask.size(); i++) {
            if (mask[i] == MASK_ZERO) {
                *coeffs[i] = 0;
            }
        }
    }
    
    void fix_mask() {
        zero_masked_coeffs();
        mask_from_coeffs();
    }
    
    virtual shared_ptr<Filter<real> > copy() {
        return make_shared<FilterSparse<real> >(filter->copy(), &mask);
    }
    
    DECLARE_COPY_APPLY(FilterSparse<real>);
    
    virtual string str(bool flatten=false) {
        shared_ptr<FilterSparse<real> > c = copy_sameclass();
        c->zero_masked_coeffs();
        return c->filter->str(flatten);
    }

    virtual void set_level(int level) {
        ASSERT(level == 0, "FilterSparse has only 1 level");
    }
    
    virtual int nlevels() {
        return 1;
    }
};

/* ------------------------------------------------------------------------
   FilterSparseMultires: Wraps any Filter in both sparse and multires
   ------------------------------------------------------------------------ */

template<class real>
class FilterSparseMultires: public FilterSparse<real> {
    public:
    shared_ptr<FilterMultires<real> > filter_multi;
    shared_ptr<Filter<real> > filter_orig;
    
    virtual ~FilterSparseMultires() { }
    
    FilterSparseMultires(const shared_ptr<Filter<real> > &filter_, vector<MaskType> *mask_=NULL, int level_=0)
      :FilterSparse<real>(shared_ptr<Filter<real> >(new FilterMultires<real>(filter_, level_)), mask_) {
        filter_multi = dynamic_pointer_cast<FilterMultires<real> >(FilterSparse<real>::filter);
        filter_orig = filter_;
        FilterSparse<real>::zero_masked_coeffs();
        FilterSparse<real>::mask_from_coeffs();
    }

    virtual shared_ptr<Filter<real> > copy() {
        return make_shared<FilterSparseMultires<real> >(filter_orig->copy(), &(FilterSparse<real>::mask), filter_multi->level);
    }
    
    //virtual DECLARE_COPY_SAMECLASS(FilterSparseMultires<real>);

    virtual void set_level(int level) {
        filter_multi->set_level(level);
        FilterSparse<real>::mask_from_coeffs();
    }

    virtual int nlevels() {
        return filter_multi->nlevels();
    }
};

/* ------------------------------------------------------------------------
   Initial parameters for FIR, IIR
   ------------------------------------------------------------------------ */

template<class real>
Array<real> fir_initial_params(int w=-1) {
    if (w < 0) { w = params.filter_w; }
    return Array<real>::random({w, w});
}

template<class real>
Array<real> fir_initial_params_h() {
    return Array<real>::random({1, params.filter_w});
}

template<class real>
Array<real> fir_initial_params_v() {
    return Array<real>::random({params.filter_w, 1});
}

template<class real>
Array<real> iir_initial_params_K(int w=-1) {
    if (w < 0) { w = params.filter_w; }
    return Array<real>::random({w});
}

template<class real>
Array<real> iir_initial_params_F(int w=-1) {
    if (w < 0) { w = params.filter_w; }
    return Array<real>::random({w});
}

#endif
