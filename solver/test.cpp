
#include "pareto.h"
#include "util.h"
#include <adept.h>

typedef ParetoPoint<int> ParetoPoint_T;
typedef Pareto<int> Pareto_T;

void print_state(const vector<ParetoPoint_T> &L, const Pareto_T &p1, const Pareto_T &p2) {
    printf("L:\n");
    for (int i = 0; i < (int) L.size(); i++) {
        printf("%f %f\n", L[i].T, L[i].err);
    }
    printf("\n");
    
    printf("p1:\n%s\n", p1.str().c_str());
    printf("p2:\n%s\n", p2.str().c_str());
}

#define ADD_AND_ASSERT(p, T, err, idx) \
{ \
    Pareto_T p0(p); \
    bool added = p.add(T, err, idx); \
    double d = p0.dist(T, err, idx); \
    ASSERT1(added ? (d < 1): (d >= 1)); \
}

void test_pareto() {
    Pareto_T p;
    ADD_AND_ASSERT(p, 1.0, 2.0, 1);  // T, err
    ADD_AND_ASSERT(p, 1.5, 2.5, 2);
    ADD_AND_ASSERT(p, 2.0, 3.0, 3);
    ADD_AND_ASSERT(p, 2.0, 1.9, 4);
    ADD_AND_ASSERT(p, 1.5, 1.8, 5);
    ADD_AND_ASSERT(p, 3.0, 2.5, 6);
    ADD_AND_ASSERT(p, 3.0, 1.7, 7);
    ADD_AND_ASSERT(p, 3.0, 1.6, 8);
    ADD_AND_ASSERT(p, 3.0, 1.6, 9);
    ADD_AND_ASSERT(p, 2.0, 3.0, 10);
    ADD_AND_ASSERT(p, 2.9, 1.6, 11);
    p.check();
    
    if (p.L.size() != 3) {
        fprintf(stderr, "Error: pareto test fail\n"); exit(1);
    }
    if (p.L[0].T != 1.0 || p.L[0].err != 2.0 || p.L[0].p != 1) {
        fprintf(stderr, "Error: pareto test 1 fail\n"); exit(1);
    }
    if (p.L[1].T != 1.5 || p.L[1].err != 1.8 || p.L[1].p != 5) {
        fprintf(stderr, "Error: pareto test 2 fail\n"); exit(1);
    }
    if (p.L[2].T != 2.9 || p.L[2].err != 1.6 || p.L[2].p != 11) {
        fprintf(stderr, "Error: pareto test 3 fail\n"); exit(1);
    }

    for (int i = 0; i < 10000; i++) {
        //printf("======== BEGIN TEST ==============\n");
        Pareto_T p1;
        int n = (rand()%100)+1;
        vector<ParetoPoint_T> L;
        double current_T = 0.5, current_E = 0.5;
        for (int j = 0; j < n; j++) {
            if (rand()%2==0) {
                current_T = rand_f();
                current_E = rand_f();
            } else if (rand()%2==0) {
                current_T += 0.1;
            } else if (rand()%2==0) {
                current_E += 0.1;
            }
            
            L.push_back(ParetoPoint_T(current_T, current_E));
            //p1.add(L[L.size()-1].T, L[L.size()-1].err, 0);
            ADD_AND_ASSERT(p1, L[L.size()-1].T, L[L.size()-1].err, 0);
        }
        //printf("========================================\n");
        Pareto_T p2;
        for (int j = n-1; j >= 0; j--) {
            //p2.add(L[j].T, L[j].err, 0);
            ADD_AND_ASSERT(p2, L[j].T, L[j].err, 0);
        }
        p1.check();
        p2.check();
        if (p1.L.size() != p2.L.size()) {
            fprintf(stderr, "p1.size() != p2.size()\n");
            print_state(L, p1, p2);
            exit(1);
        }
        for (int j = 0; j < (int) p1.L.size(); j++) {
            if (!(p1.L[j] == p2.L[j])) {
                fprintf(stderr, "Error: pareto test 2 fail\n");
                print_state(L, p1, p2);
                exit(1);
            }
        }
    }
    
    printf("Pareto:        OK\n");
}

int main () {
    test_pareto();
    
    return 0;
}

