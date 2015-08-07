
#ifndef _util_h
#define _util_h

#include <string>
#include <memory>
#include <vector>
#include <cassert>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

using std::vector;
using std::string;
using std::make_shared;

/* Things from CImg */
inline double log2(const double x) {
  static const double base = std::log(2.0);
  return std::log(x)/base;
}


/* Assert that never compiles away */
#define ASSERT1(x) if (!(x)) { fprintf(stderr, "Assert failed at %d\n", __LINE__); assert(false); exit(1); }
#define ASSERT2(x, msg) if (!(x)) { fprintf(stderr, "%s\n", msg); fprintf(stderr, "Assert failed at %d\n", __LINE__); assert(false); exit(1); }

string str_replace(string s, const string &toReplace, const string &replaceWith);

bool file_exists(const char *filename);

vector<string> str_split(const string &s, char c);

string lstrip(string s);
string rstrip(string s);
string strip(string s);
int count(string s, char sub);

string read_binary_file(const string &filename);

string exec_output(const char *cmd);

class ProbabilitySampler { public:
    vector<double> C;
    ProbabilitySampler(const vector<double> &prob);
    int sample();
};

/* Trapezoid rule integration */
double integrate_trapezoid(const vector<double> &x, const vector<double> &y);

/* Integration using boxes that contain the maximum y value */
double integrate_box_max(const vector<double> &x, const vector<double> &y);

void logf(const char *fmt, ...);

double product(const vector<int> &L);

vector<int> multiply(const vector<int> &L, double factor);

double sinc(double x);

template<class real>
real max(const vector<real> &L) {
    ASSERT2(L.size(), "expected nonzero length for max()");
    real ans(L[0]);
    for (int i = 1; i < (int) L.size(); i++) {
        if (L[i] > ans) {
            ans = L[i];
        }
    }
    return ans;
}

#endif

