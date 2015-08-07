
#include "iadouble.h"

using namespace adept;

namespace std {
    using ::cos;
    using ::sin;
    using ::exp;
}

#include <boost/numeric/interval.hpp>

adouble square(adouble x) {
    return x * x;
}

adouble hull(adouble x, adouble y) {
    ASSERT2(false, "hull() should not be called with adouble");
    return 0.0;
}

double hull(double x, double y) {
    ASSERT2(false, "hull() should not be called with double");
    return 0.0;
}

