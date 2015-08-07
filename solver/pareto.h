
/* Pareto frontier. Tracks a monotonically non-increasing function E = f(T) made by point tuples (T, E).  */

#ifndef _pareto_h
#define _pareto_h

#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include "array.h"
#include "timer.h"
#include "util.h"
#include "filter.h"
#include "params.h"

using std::string;
using std::vector;
using std::map;
using std::shared_ptr;

template<class C>
class ParetoPoint { public:
    double T, err;
    C p;
    ParetoPoint(double T_=0, double err_=0, C p_=C())
     :T(T_), err(err_), p(p_) { }
};

template<class C>
bool operator <  (const ParetoPoint<C> &a, const ParetoPoint<C> &b) {
    if (a.T < b.T) {
        return true;
    } else if (a.T > b.T) {
        return false;
    }
    return a.err < b.err;
}

template<class C>
bool operator == (const ParetoPoint<C> &a, const ParetoPoint<C> &b) {
    return a.T == b.T && a.err == b.err;
}
    
template<class C>
class Pareto { public:
    vector<ParetoPoint<C> > L;
    map<int, int> points_per;
    int max_points_per;
    double Tgrid;
    
    Pareto(int max_points_per_=1000*1000, double Tgrid_=1.0)
      :max_points_per(max_points_per_), Tgrid(Tgrid_) { }

    void clear() {
        L.clear();
        points_per.clear();
    }
    
    void check() {
        for (int i = 0; i < (int) (L.size()-1); i++) {
            ParetoPoint<C> p1(L[i]);
            ParetoPoint<C> p2(L[i+1]);
            if (p1.T >= p2.T) {
                fprintf(stderr, "p1.T (%17e) >= p2.T (%17e) at i=%d, p1.err=%.17e, p2.err=%.17e\n", p1.T, p2.T, i, p1.err, p2.err);
                ASSERT(0, "Pareto::check failed\n");
            }
            if (p2.err >= p1.err) {
                fprintf(stderr, "p2.err (%17e) >= p1.err (%17e) at i=%d, p1.T=%.17e, p2.T=%.17e\n", p2.err, p1.err, i, p1.T, p2.T);
                ASSERT(0, "Pareto::check failed\n");
            }
        }
    }

    double dist(double T, double err, C p) {   /* Returns distance ratio >= 1 if outside Pareto Frontier, else < 1. */
        ParetoPoint<C> pt(T, err, p);
        int i = std::upper_bound(L.begin(), L.end(), pt) - L.begin();
        //fprintf(stderr, "Pareto::add(%.17e, %.17e), i=%d, L[i-1]=(%.17e, %.17e)\n", T, err, i, i-1 >= 0 ? L[i-1].T: 0.0, i-1>=0 ? L[i-1].err: 0.0);
        double ratio = 0.0;
        if (i - 1 >= 0) {
            ratio = MIN(err/L[i-1].err, T/L[i-1].T);
        }
        
        return ratio;
    }
    
    bool add(double T, double err, C p)  {
        ParetoPoint<C> pt(T, err, p);
        int i = std::upper_bound(L.begin(), L.end(), pt) - L.begin();
        //fprintf(stderr, "add(), current Pareto:\n%s\n", str().c_str());
        //fprintf(stderr, "Pareto::add(%.17e, %.17e), i=%d, L[i-1]=(%.17e, %.17e)\n\n\n", T, err, i, i-1 >= 0 ? L[i-1].T: 0.0, i-1>=0 ? L[i-1].err: 0.0);
        if (i - 1 >= 0 && T >= L[i-1].T && err >= L[i-1].err) {
            return false;
        }
        int gridv = int(T/Tgrid);
        if (points_per[gridv] >= max_points_per) {
            return false;
        }
        
        points_per[gridv]++;
        L.insert(L.begin()+i, pt);
        
        //printf("Before restore invariant:\n%s\n\n", str().c_str());
        
        /* Restore invariant */
        auto iter = L.begin();
        ParetoPoint<C> prev;
        bool prev_exists = false;
        while (iter != L.end()) {
            auto current = *iter;
            if (prev_exists) {
                if (current.err >= prev.err) {
                    int gridv_p = int(current.T/Tgrid);
                    points_per[gridv_p]--;
                    ASSERT(points_per[gridv_p] >= 0, "points_per < 0");
                    iter = L.erase(iter);
                    continue;
                }
            }
            prev_exists = true;
            prev = current;
            
            ++iter;
        }
        
        //printf("After restore invariant:\n%s\n\n", str().c_str());
        
        //check();
        return true;
    }    
        
    string str() const {
        string ans("");
        char buf[256];
        for (int i = 0; i < (int) L.size(); i++) {
            sprintf(buf, "%f %f\n", L[i].T, L[i].err);
            ans += buf;
        }
        return ans;
    }
    
    void join(const Pareto<C> &other) {
        for (int i = 0; i < (int) other.L.size(); i++) {
            add(other.L[i].T, other.L[i].err, other.L[i].p);
        }
    }
    
    bool check_min_points() {
        int nsuccess = 0;
        for (int i = 0; i < (int) L.size(); i++) {
            double err = L[L.size()-1-i].err;
            if (err < params.pareto_min_points_error) {
                nsuccess++;
            }
        }
        return (nsuccess >= params.pareto_min_points);
    }
    
    double max_time(double max_E=1.0) {
        ASSERT(L.size(), "max_time() requires non-empty pareto");
        if (!check_min_points()) { return 1e6; }
        return L[L.size()-1].T;
    }
    
    double avg_time(double max_E=1.0) {
        ASSERT(L.size(), "avg_time() requires non-empty pareto");
        if (!check_min_points()) { return 1e6; }
        
        vector<double> T;
        vector<double> E;
        E.push_back(0.0);
        double Tval = L[L.size()-1].T;
        T.push_back(Tval);
        for (int i = 0; i < (int) L.size(); i++) {
            double err = L[L.size()-1-i].err;
            if (err < max_E) {
                Tval = L[L.size()-1-i].T;
                E.push_back(err);
                T.push_back(Tval);
            }
        }
        
        E.push_back(max_E);
        T.push_back(Tval);
        return integrate_box_max(E, T) / max_E;
//        return integrate_trapezoid(E, T) / max_E;
    }
};

template<class real>
class ParetoSparse { public:
    typedef Pareto<shared_ptr<FilterSparse<real> > > type;
};

#define PARETO_TYPE(real) typename ParetoSparse<real>::type

template<class real>
class ParetoTraceWriter { public:
    string filename;
    bool init;
    double T0;
    PARETO_TYPE(real) base;
    vector<double> Tavg_L;
    ParetoTraceWriter(string filename_, double T0_, const Array<real> &in, Array<real> &out, const Array<double> &target) :filename(filename_), init(false), T0(T0_) {
        if (params.merge_base) {
            get_base_solution(in, out, target, &base);
        }
    }
    
    void write(PARETO_TYPE(real) *pareto0, int total_of_calls=0) {
        //double T_start = wall_time();
        PARETO_TYPE(real) *pareto = pareto0;
        PARETO_TYPE(real) pcopy;
        if (params.merge_base) {
            pcopy = *pareto0;
            pcopy.join(base);
            pareto = &pcopy;
        }
        //printf("merge_base=%d, pareto size: %d, base size: %d, pcopy size: %d\n", params.merge_base, (int) pareto0->L.size(), (int) base.L.size(), (int) pcopy.L.size());
        FILE *f = fopen(filename.c_str(), init ? "at": "wt");
        if (!f) { fprintf(stderr, "Error opening '%s'\n", filename.c_str()); exit(1); }
        if (!init) { fprintf(f, "["); }
        else { fprintf(f, ",\n"); }
        double avgT = pareto->avg_time();
        Tavg_L.push_back(avgT);
        fprintf(f, "{\"wall_time\": %f, \"of_calls\": %d, \"avg_time\": %f, \"avg_time10\": %f, \"avg_time20\": %f, \"pareto\": [", wall_time()-T0, total_of_calls, avgT, pareto->avg_time(0.1), pareto->avg_time(0.2));
        for (int i = 0; i < (int) pareto0->L.size(); i++) {
            fprintf(f, "[%f, %f]", pareto0->L[i].T, pareto0->L[i].err);
            if (i < ((int) pareto0->L.size())-1) { fprintf(f, ",\n"); }
        }
        fprintf(f, "]}");
        fclose(f);
        init = true;
        //double T_end = wall_time();
        //printf("write pareto trace time: %f secs\n", T_end-T_start);
    }
    ~ParetoTraceWriter() {
//        printf("destructor\n");
        if (init) {
//            printf("destructor and init, filename: %s\n", filename.c_str());
            FILE *f = fopen(filename.c_str(), "at");
            fprintf(f, "]");
            fclose(f);
        }
    }
};

template<class real>
void write_pareto_trace(PARETO_TYPE(real) *pareto, string filename, double T0, int total_of_calls, const Array<real> &in, Array<real> &out, const Array<double> &target, vector<double> *Tavg_L_out) {
    bool free = false;
    static map<string, shared_ptr<ParetoTraceWriter<real> > > d;
    if (free) {
//        printf("calling free\n");
        d.erase(filename);
        return;
    }
    if (!d.count(filename)) {
        d[filename] = make_shared<ParetoTraceWriter<real> >(filename, T0, in, out, target);
    }
    d[filename]->write(pareto, total_of_calls);
    if (Tavg_L_out) {
        *Tavg_L_out = d[filename]->Tavg_L;
    }
}

#endif
