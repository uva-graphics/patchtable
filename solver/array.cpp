#include "array.h"

//bool array_debug = false;

#define COPY_ARRAY(source) \
void copy_array(const Array<source> &src, Array<source> &dest) { \
    dest = src; \
}

COPY_ARRAY(double);
COPY_ARRAY(float);
COPY_ARRAY(uint8_t);
COPY_ARRAY(uint16_t);
COPY_ARRAY(uint32_t);
COPY_ARRAY(int8_t);
COPY_ARRAY(int16_t);
COPY_ARRAY(int32_t);

#define real float

vector<real> gaussian_kernel(real sigma) {
    sigma = std::abs(sigma);
    
    real q;
    
    if (sigma >= 2.5) {
        q = 0.98711 * sigma - 0.96330;
    } else {
        q = 3.97156 - 4.14554 * sqrt(1 - 0.26891*sigma);
        q = MAX(q, 0.0);
    }
    
    real b0 = 1.57825 + 2.44413 * q + 1.4281 * q*q + 0.422205 * q*q*q;
    real b1 = 2.44413 * q + 2.85619 * q*q + 1.26661 * q*q*q;
    real b2 = -(1.4281 * q*q + 1.26661 * q*q*q);
    real b3 = 0.422205 * q*q*q;
    
    real B = 1 - (b1 + b2 + b3) / b0;
    
    vector<real> tempV;
    tempV.push_back(B);
    tempV.push_back(b1/b0);
    tempV.push_back(b2/b0);
    tempV.push_back(b3/b0);
    tempV.push_back(0.0);


    return tempV;
}

string get_additional_filename(const char *filename, int n) {
    string filename_p(filename);
    filename_p = filename_p.substr(0, filename_p.size()-4);
    char buf[256];
    sprintf(buf, "%d.pfm", n);
    filename_p += buf;
    return filename_p;
}

bool string_ends_with(const char *a0, const char *b0) {
    string a(a0);
    string b(b0);
    if (a.length() >= b.length()) {
        return a.compare(a.length()-b.length(), b.length(), b) == 0;
    }
    return false;
}
