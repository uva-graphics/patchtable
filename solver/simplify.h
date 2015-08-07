#ifndef _simplify_h
#define _simplify_h

#include "filter.h"
#include "util.h"
#include <set>
#include <map>

using std::map;
using std::set;

#define CLEAN_SYMBOLS 1

template<class real>
class SimplifyRecurse { public:
    set<int> seen;
    map<int, string> idx_str;
    int max_idx;
    
    SimplifyRecurse() {
        max_idx = 0;
    }
    
    string index_to_string(int i) {
        if (idx_str.count(i)) {
            return idx_str[i];
        }
        idx_str[i] = string("p[") + to_string(max_idx) + string("]");
        max_idx++;
        return idx_str[i];
    }
    
    string simplify_recurse(const shared_ptr<FilterDAG<real> > &f, int current) {
        if (current == INPUT_SOURCE_IMAGE) {
            return "Input";
        }
        ASSERT(current >= 0 && current < (int) f->L.size(), "DAG node out of bounds");
        shared_ptr<Filter<real> > fcurrent = f->L[current].f;
        string f_name = string(typeid(*fcurrent).name());
        
        if (seen.count(current)) {
            return index_to_string(current);
        }
        vector<int> inputs = f->L[current].inputs;
        if (typeid(*fcurrent) == typeid(FilterAdd<real>)) {
            //printf("Encountered FilterAdd\n");
            int i = 0;
            while (i < (int) inputs.size()) {
                ASSERT(i >= 0 && i < (int) inputs.size(), "index i out of bounds");
                bool changed = false;
                if (inputs[i] != INPUT_SOURCE_IMAGE) {
                    ASSERT(inputs[i] >= 0 && inputs[i] < (int) f->L.size(), "input out of bounds");
                    NodeDAG<real> input = f->L[inputs[i]];
                    if (typeid(*input.f) == typeid(FilterAdd<real>)) {
                        //printf("  Found FilterAdd inside FilterAdd\n");
                        inputs.erase(inputs.begin()+i);
                        inputs.insert(inputs.begin()+i, input.inputs.begin(), input.inputs.end());
                        changed = true;
                    }
                }
                if (!changed) {
                    i++;
                }
            }
            sort(inputs.begin(), inputs.end());
        }
        seen.insert(current);
        string ans = string(f_name) + "[" + index_to_string(current) + "](";
        for (int i = 0; i < (int) inputs.size(); i++) {
            ans += simplify_recurse(f, inputs[i]);
            if (i < int(inputs.size())-1) { ans += ", "; }
        }
        ans += ")";
        return ans;
    }
};

template<class real>
string simplify(const shared_ptr<FilterDAG<real> > &f) {
    f->check();
//    printf("Simplify:\n%s\n\n", f->str().c_str());
    ASSERT(f->L.size(), "simplify on length 0 DAG");
    //if (!f->L.size()) { return "(Empty DAG)"; }
    string ans = SimplifyRecurse<real>().simplify_recurse(f, f->L.size()-1);

#if CLEAN_SYMBOLS
    ans = str_replace(ans, "IN5adept5aRealEE", "");
    ans = str_replace(ans, "25Filter", "Filter");
    ans = str_replace(ans, "24Filter", "Filter");
    ans = str_replace(ans, "16Filter", "Filter");
    ans = str_replace(ans, "15Filter", "Filter");
    ans = str_replace(ans, "14Filter", "Filter");
    ans = str_replace(ans, "9Filter", "Filter");
#endif
    return ans;
}

template<class real>
shared_ptr<FilterDAG<real> > program1() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 0}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program2() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 1}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program3() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {1, 0}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program4() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {1, 0}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program5() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {1, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {2, 0}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program6() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 1}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 2}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program7() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {0, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {2, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {3, 1}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
shared_ptr<FilterDAG<real> > program8() {
    vector<NodeDAG<real> > L;
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterFIR<real>(Array<real>())), {INPUT_SOURCE_IMAGE}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {1, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {2, 0}));
    L.push_back(NodeDAG<real>(shared_ptr<Filter<real> >(new FilterAdd<real>()), {3, 0}));
    return shared_ptr<FilterDAG<real> >(new FilterDAG<real>(L, false));
}

template<class real>
void test_simplify() {
    string sp1 = simplify(program1<real>());
    string sp2 = simplify(program2<real>());
    string sp3 = simplify(program3<real>());
    string sp4 = simplify(program4<real>());
    string sp5 = simplify(program5<real>());
    string sp6 = simplify(program6<real>());
    string sp7 = simplify(program7<real>());
    string sp8 = simplify(program8<real>());
    //printf("sp3:\n%s\n\n", sp3.c_str());
    //printf("sp4:\n%s\n\n", sp4.c_str());
    
    ASSERT1(sp1 != sp2);
    ASSERT1(sp3 == sp2);
    ASSERT1(sp4 == sp3);
    ASSERT1(sp5 != sp4);
    ASSERT1(sp5 != sp3);
    ASSERT1(sp5 != sp2);
    ASSERT1(sp5 != sp1);
    ASSERT1(sp6 == sp5);
    ASSERT1(sp7 != sp1);
    ASSERT1(sp7 != sp2);
    ASSERT1(sp7 != sp3);
    ASSERT1(sp7 != sp4);
    ASSERT1(sp7 != sp5);
    ASSERT1(sp7 != sp6);
    ASSERT1(sp8 == sp7);

    printf("simplify:   OK\n");
}

#endif
