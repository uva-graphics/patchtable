/* Adept double using interval arithmetic (via boost::numeric::interval) */

#ifndef _iadouble_h
#define _iadouble_h

#include <adept.h>
#include "util.h"

using namespace adept;

namespace std {
    using ::cos;
    using ::sin;
    using ::exp;
}

#include <boost/numeric/interval.hpp>

using namespace boost::numeric::interval_lib;
using boost::numeric::interval;

typedef interval<adouble, policies<save_state_nothing<rounded_transc_exact<adouble> >, checking_base<adouble> > > iadouble;

adouble square(adouble x);
adouble hull(adouble x, adouble y);
double hull(double x, double y);

#endif
