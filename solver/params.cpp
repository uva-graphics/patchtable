
#include "params.h"
#include "array.h"

Params params;

vector<double> *Params::feature_coeff_list() {
    if (!feature_coeff_list_internal.size()) {
        Array<double> M = load_matrix<double>(feature_coeff);
        feature_coeff_list_internal.resize(M.nelems);
        for (int i = 0; i < M.nelems; i++) {
            feature_coeff_list_internal[i] = M.data[i];
        }
    }
    return &feature_coeff_list_internal;
}

Array<double> *Params::weights_array() {
    if (!weights.size()) { return NULL; }
    if (weights_array_internal.dimensions() < 2) {
        weights_array_internal = load_matrix<double>(weights);
        if (params.target_w > 0) {
            weights_array_internal = zero_pad_center(weights_array_internal, params.target_w, params.target_w);
        }
    }
    return &weights_array_internal;
}
