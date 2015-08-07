
#include "filter.h"

bool log_dag = false;
double T_fir = 0.0;

int downsample_dim(int w) {
    w = w/2;
    if (w % 2 == 0) { w++; }
    return w;
}

Array<int> downsample_box(const Array<int> &box) {
    ASSERT(box.dimensions() == 2, "expected 2D box");
    int h0 = box.height(), w0 = box.width();
    int h = downsample_dim(h0), w = downsample_dim(w0);
    if (h <= 1 || w <= 1) {
        return box;
    }
    return resample_nearest(box, {h, w}, true);
}

vector<int> flatten_boxes(const boxes_t &L) {
    vector<int> ans;
    for (int i = 0; i < (int) L.size(); i++) {
        shared_ptr<Array<int> > box(L[i]);
        for (int j = 0; j < box->nelems; j++) {
            ans.push_back(box->data[j]);
        }
    }
    return ans;
}

boxes_t copy_boxes(const boxes_t &b) {
    boxes_t ans;
    for (int i = 0; i < (int) b.size(); i++) {
        ans.push_back(make_shared<Array<int> >(*b[i]));
    }
    return ans;
}

void test_multiscale() {
    Array<int> A;
    int h = 15; //11;
    int w = 15; //10;
    A.resize({h, w});
    for (int i = 0; i < A.nelems; i++) {
        A.data[i] = i;
    }
    ASSERT(equals(A, resample_nearest(A, A.sizes, true)), "resample(A.sizes) != A");
    ASSERT(equals(A, resample_nearest(A, A.sizes, false)), "resample(A.sizes) != A");
    for (int j = 0; j < 7; j++) {
        printf("%s\n\n", A.str().c_str());
        printf("%s\n\n", resample_nearest(A, {h, w}, true).str().c_str());
        printf("---------------------------------------------\n");
        A = downsample_box(A);
    }
}

MaskType int_to_masktype(int x) {
    if (x == 0) {
        return MASK_ZERO;
    } else if (x == 1) {
        return MASK_VAR;
    } else {
        ASSERT(false, "int_to_masktype() received invalid argument");
        fprintf(stderr, "int_to_masktype() failed\n");
        exit(1);
    }
}