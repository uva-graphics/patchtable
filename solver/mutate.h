
#ifndef _mutate_h
#define _mutate_h

#include "filter.h"
#include "simplify.h"
#include <algorithm>

using std::pair;

#define MAX_MUTATE_ITERS         1000
#define POINT_MUTATE_DOWNSAMPLE  1

/* Insert node to DAG, before existing node at position i, making old node read new node as input. */
template<class real>
void insert_dag(shared_ptr<FilterDAG<real> > f, int i, shared_ptr<Filter<real> > insert_f, const vector<int> &insert_inputs, bool correct_all=false) {
    ASSERT(i >= 0 && i <= (int) f->L.size(), "insert_dag: index out of bounds");
    f->L.insert(f->L.begin()+i, NodeDAG<real>(insert_f, insert_inputs));
    
    /* Fix indices in DAG. */
    for (int j = i+1; j < (int) f->L.size(); j++) {
        for (int k = 0; k < (int) f->L[j].inputs.size(); k++) {
            if (j == i + 1 && (k == 0 || correct_all)) {         /* Connect old node to read new node as input */
                f->L[j].inputs[k] = i;
                continue;
            }
            if (f->L[j].inputs[k] >= i) {
                ASSERT(f->L[j].inputs[k] != INPUT_SOURCE_IMAGE, "input should not be source image");
                f->L[j].inputs[k]++;
                ASSERT(f->L[j].inputs[k] <= (int) f->L.size(), "adjusted input out of range");
            }
        }
    }
}

template<class real>
void delete_dag(shared_ptr<FilterDAG<real> > f, int i) {
    ASSERT(i >= 0 && i < (int) f->L.size(), "delete_dag: index out of bounds");
    ASSERT(f->L[i].inputs.size(), "expected delete_dag node to have an input");
    int forward_idx = f->L[i].inputs[0];
    f->L.erase(f->L.begin()+i);
    
    /* Fix DAG indices, and also forward inputs >= i to the input of the removed node. */
    for (int j = i; j < (int) f->L.size(); j++) {
        for (int k = 0; k < (int) f->L[j].inputs.size(); k++) {
            if (f->L[j].inputs[k] == i) {
                f->L[j].inputs[k] = forward_idx;
            } else if (f->L[j].inputs[k] > i) {
                f->L[j].inputs[k]--;
            }
        }
    }
}

template<class real>
bool add_node(shared_ptr<FilterDAG<real> > f, set<string> &seen, vector<shared_ptr<FilterDAG<real> > > &ans) {
    string s = simplify(f);
    if (params.verbose >= 2) {
        printf("Filter is:\n%s\n\n", f->str().c_str());
        printf("Simplify gives %d: %s\n", (int) f->L.size(), s.c_str());
        printf("Seen is: %d\n", (int) seen.count(s));
    }
    if (!seen.count(s)) {
        seen.insert(s);
        ans.push_back(f);
        return true;
    }
    return false;
}

vector<int> random_order(int n, bool is_random=true);
vector<int> random_order(int lo, int hi, bool is_random=true);

#define RANDOM_ORDER1(n) random_order((n), is_random)
#define RANDOM_ORDER2(lo, hi) random_order((lo), (hi), is_random)

#define ADD_NODE(msg) \
    if (f->check_sizes(input_size)) { \
        add_node(f, seen, ans); \
        if (!all && ans.size()) { return ans; } \
    }

// if (string(msg).size()) { printf("%s\n", msg); } 

template<class real>
vector<int> caller_distance(shared_ptr<FilterDAG<real> > f0, int jnode) {
    vector<int> dist(int(f0->L.size()), -1);          /* Distance of callers (and recursive callers of those) from node j */
    ASSERT((unsigned) jnode < f0->L.size(), "jnode out of bounds in caller_distance");
    dist[jnode] = 0;
    for (int i = jnode+1; i < (int) f0->L.size(); i++) {
        for (int k = 0; k < (int) f0->L[i].inputs.size(); k++) {
            int input = f0->L[i].inputs[k];
            if (input >= 0 && (input == jnode || dist[input] > 0)) {
                dist[i] = MAX(dist[input], 0) + 1;
            }
        }
    }
    return dist;
}

template<class real>
int get_upsample_pair(shared_ptr<FilterDAG<real> > f0, int jnode) {
    vector<int> dist = caller_distance(f0, jnode);
            
    int least_dist_idx = -1;                         /* UpsamplePostfilter node index with least distance */
    int least_dist = 1000*1000;
    for (int i = jnode+1; i < (int) f0->L.size(); i++) {
        if (typeid(*f0->L[i].f) == typeid(FilterUpsamplePostfilter<real>) && dist[i] >= 0 && dist[i] < least_dist) {
            least_dist = dist[i];
            least_dist_idx = i;
        }
    }

    return least_dist_idx;
}

/* Get one or all mutations that remove one filter from the DAG (or two filters for paired filters). */
template<class real>
vector<shared_ptr<FilterDAG<real> > > delete_mutate(shared_ptr<FilterDAG<real> > f0, set<string> &seen, bool all=false) {
    vector<shared_ptr<FilterDAG<real> > > ans;
    vector<double> input_size;
    bool is_random = !all;

    if (f0->L.size() == 0) { return {}; }
    
//    int jchosen = rand()%int(f0->L.size());
    for (int jnode: RANDOM_ORDER1(int(f0->L.size()))) {      /* Node to delete */
//        if (!all && jnode != jchosen) { continue; }
        
        auto f = f0->copy_sameclass();
        if (typeid(*f0->L[jnode].f) == typeid(FilterDownsamplePrefilter<real>)) {
            int least_dist_idx = get_upsample_pair(f0, jnode);
            if (least_dist_idx < 0) {       /* No upsample caller (probably some dead code should have been eliminated) */
                continue;
            }

            ASSERT(least_dist_idx >= 0, "could not find upsample caller");
            delete_dag(f, jnode);
            ASSERT((unsigned) (least_dist_idx-1) < (unsigned) f->L.size(), "least_dist_idx out of bounds after removing jnode");
            delete_dag(f, least_dist_idx-1);
            if (f->L.size()) {
                ADD_NODE("delete_mutate success for DownsamplePrefilter");
            } else {
                return {};
            }
        } else {
            delete_dag(f, jnode);
            if (f->L.size()) {
                ADD_NODE("delete_mutate success for basic");
            } else {
                return {};
            }
        }
    }
    return ans;
}

template<class real>
void get_basic_filters(vector<shared_ptr<Filter<real> > > &fadd, vector<int> &nargs) {
    // FilterFIR
    fadd.push_back(shared_ptr<Filter<real> >(new FilterFIR<real>(fir_initial_params<real>())));
    nargs.push_back(1);

    if (params.vh_mode) {
        // FilterFIRH
        fadd.push_back(shared_ptr<Filter<real> >(new FilterFIRH<real>(fir_initial_params_h<real>())));
        nargs.push_back(1);

        // FilterFIRV
        fadd.push_back(shared_ptr<Filter<real> >(new FilterFIRV<real>(fir_initial_params_v<real>())));
        nargs.push_back(1);
    }

    // FilterIIR (+x direction)
    fadd.push_back(shared_ptr<Filter<real> >(new FilterIIR<real>(iir_initial_params_K<real>(), iir_initial_params_F<real>())));
    nargs.push_back(1);

    // FilterIIR (-x direction)
    fadd.push_back(shared_ptr<Filter<real> >(new FilterIIRMinusX<real>(iir_initial_params_K<real>(), iir_initial_params_F<real>())));
    nargs.push_back(1);

    // FilterIIR (+y direction)
    fadd.push_back(shared_ptr<Filter<real> >(new FilterIIRPlusY<real>(iir_initial_params_K<real>(), iir_initial_params_F<real>())));
    nargs.push_back(1);

    // FilterIIR (-y direction)
    fadd.push_back(shared_ptr<Filter<real> >(new FilterIIRMinusY<real>(iir_initial_params_K<real>(), iir_initial_params_F<real>())));
    nargs.push_back(1);
}

/* Insert subgraph b into DAG a, just prior to idx_a. */
template<class real>
int insert_subgraph(shared_ptr<FilterDAG<real> > a, shared_ptr<FilterDAG<real> > b, int idx_a, int orig_input) {
    b->check();
    ASSERT((unsigned) idx_a <= a->L.size(), "idx_a out of bounds");
    a->L.insert(a->L.begin()+idx_a, b->L.begin(), b->L.end());
    int insert_end = idx_a + int(b->L.size());
    ASSERT(insert_end <= (int) a->L.size(), "insert_end out of bounds");
    
    /* Fix indices in DAG. */
    for (int j = idx_a; j < (int) a->L.size(); j++) {
        for (int k = 0; k < (int) a->L[j].inputs.size(); k++) {
            if (a->L[j].inputs[k] < 0) {
                a->L[j].inputs[k] = orig_input;     /* Point -1 to the input to the node we inserted before. */
                continue;
            }
            if (j == insert_end) {         /* Connect node inserted before to read output of b subgraph as input */
                ASSERT(a->L[j].inputs.size() == 1, "expected 1 input");
                a->L[j].inputs[k] = insert_end-1;
                continue;
            }
            int delta = 0;
            
            if (j > insert_end && a->L[j].inputs[k] >= idx_a) {
                delta = int(b->L.size());
            } else if (j < insert_end) {
                delta = idx_a;
            }
            //if (a->L[j].inputs[k] >= idx_a) { //insert_end) {
            a->L[j].inputs[k] += delta; //(j > insert_end) ? int(b->L.size()): idx_a;
            ASSERT(a->L[j].inputs[k] < (int) a->L.size(), "adjusted input out of range");
            ASSERT(a->L[j].inputs[k] < j, "adjusted input exceeds current node");
            //}
        }
    }
    return insert_end;
}

/* Mutation similar to Convolution Pyramids paper. Replace 'H' with:
 
        ----------------------------------> H -----------------------------> + --->  Out
        \                                                                 /
         \--> (Prefilter -> Downsample) -> F -> (Upsample -> Postfilter) /

   Note that FilterDownsamplePrefilter and FilterUpsamplePostfilter are used to ensure shift-invariance and keep
   the upsample/downsample from being split from the pre- or post-filters.
   Here H is the last subgraph* contained in a downsample-upsample block, or if no downsample-upsample
   block exists, the entire graph. Also F is an FIR or one of the IIR filters chosen at random.
   
   * The subgraph needs to be disconnected from the rest of the graph. */
   
template<class real>
vector<shared_ptr<FilterDAG<real> > > add_mutate_band(shared_ptr<FilterDAG<real> > f0, set<string> &seen, bool force_fir=false) {
    if (params.verbose >= 2) {
        printf("add_mutate_band\n");
    }
    if (f0->L.size() == 0) { return {}; }
    bool all = false;
    vector<double> input_size;
    vector<shared_ptr<FilterDAG<real> > > ans;
    int downsample_idx = -2;
    int upsample_idx = -2;
    for (int i = 0; i < (int) f0->L.size(); i++) {
        if (typeid(*f0->L[i].f) == typeid(FilterDownsamplePrefilter<real>)) {
            int up_idx = get_upsample_pair(f0, i);
            if (params.verbose >= 2) {
                printf("get_upsample_pair(%d) = %d\n", i, up_idx);
            }
            if (up_idx >= 0) {
                downsample_idx = i;
                upsample_idx = up_idx;
            }
        }
    }
    
    if (upsample_idx < 0) {
        downsample_idx = -1;
        upsample_idx = f0->L.size();
    }

    int upsample_prev = f0->L.size()-1;                 // Last node of H
    if ((unsigned) upsample_idx < f0->L.size()) {
        ASSERT(f0->L[upsample_idx].inputs.size() == 1, "inputs not size 1");
        upsample_prev = f0->L[upsample_idx].inputs[0];
    }

    vector<shared_ptr<Filter<real> > > fadd_L1;
    vector<shared_ptr<Filter<real> > > fadd_L2;
    vector<shared_ptr<Filter<real> > > fadd_L3;
    vector<int> nargs;
    get_basic_filters(fadd_L1, nargs);
    get_basic_filters(fadd_L2, nargs);
    get_basic_filters(fadd_L3, nargs);
    
    int fadd_i = rand()%fadd_L1.size();
    if (force_fir) { fadd_i = 0; }
    
    ASSERT(downsample_idx < upsample_idx, "expected downsample_idx < upsample_idx");
    auto node_a = shared_ptr<Filter<real> >(new FilterDownsamplePrefilter<real>());
    auto node_b = fadd_L2[fadd_i]->copy();
    auto node_c = shared_ptr<Filter<real> >(new FilterUpsamplePostfilter<real>());
    auto node_d = shared_ptr<Filter<real> >(new FilterAdd<real>());
    auto node_L = vector<shared_ptr<Filter<real> > >({node_a, node_b, node_c, node_d});
    auto G = make_shared<FilterDAG<real> >(node_L);
    
    auto f = f0->copy_sameclass();

    int index_last = insert_subgraph(f, G, upsample_idx, downsample_idx) - 1;
    ASSERT((unsigned) index_last < f->L.size(), "index_last out of bounds");
    
    /* Add H as an input to + */
    ASSERT(f->L[index_last].inputs.size() == 1, "expected index_last to have 1 input");
    ASSERT(typeid(*f->L[index_last].f) == typeid(FilterAdd<real>), "expected index_last to be FilterAdd");
    f->L[index_last].inputs.push_back(upsample_prev);
    
    if (params.verbose >= 2) {
        printf("add_mutate_band:\n");
        printf("downsample_idx: %d, upsample_idx: %d, upsample_prev: %d\n", downsample_idx, upsample_idx, upsample_prev);
        printf("Original:\n%s\n\nResult:\n%s\n", f0->str().c_str(), f->str().c_str());
    }
    f->check();
    
    ADD_NODE("add_mutate_band success");
    
    return {};
}

/* Replacement mutation, returns list of size zero or one. */
template<class real>
vector<shared_ptr<FilterDAG<real> > > replace_mutate(shared_ptr<FilterDAG<real> > f0, set<string> &seen) {
    bool all = false;
    vector<double> input_size;
    vector<shared_ptr<Filter<real> > > fadd;
    vector<int> nargs;
    get_basic_filters(fadd, nargs);

    vector<shared_ptr<FilterDAG<real> > > ans;
    for (int i: random_order(int(f0->L.size()))) { //for (int i = 0; i < (int) f0->size(); i++) {
        bool do_mutate = false;
        for (int j = 0; j < (int) fadd.size(); j++) {
            if (typeid(*f0->L[i].f) == typeid(*fadd[j])) {
                ASSERT(f0->L[i].inputs.size() == 1, "expected 1 input");
                do_mutate = true;
                break;
            }
        }
        if (do_mutate) {
            auto f = f0->copy_sameclass();
            f->L[i].f = fadd[rand()%fadd.size()]->copy();
            ADD_NODE("replace_mutate success");
        }
    }
    return ans;
}

/* Get one or all mutations that add one filter to the DAG (or two filters for paired filters). */
template<class real>
vector<shared_ptr<FilterDAG<real> > > add_mutate(shared_ptr<FilterDAG<real> > f0, set<string> &seen, bool all=false) {
    // Unary:  IIR/FIR
    // Binary: FilterAdd
    // Paired: DownsamplePrefilter/UpsamplePostfilter, Transpose/Transpose, FlipH/FlipH

//    ASSERT(target.dimensions() == 2, "target dimensions not 2");
    vector<shared_ptr<FilterDAG<real> > > ans;
    
    vector<shared_ptr<Filter<real> > > fadd;
    vector<int> nargs;
    get_basic_filters(fadd, nargs);

    // FilterAdd
    fadd.push_back(shared_ptr<Filter<real> >(new FilterAdd<real>()));
    nargs.push_back(2);

#if POINT_MUTATE_DOWNSAMPLE
    // DownsamplePrefilter/UpsamplePostfilter
    if (f0->L.size() && params.ga_downsample) {
        fadd.push_back(shared_ptr<Filter<real> >(new FilterDownsamplePrefilter<real>()));
        nargs.push_back(1);
    }
#endif

    vector<double> input_size;
    bool is_random = !all;
    
    ASSERT(fadd.size() == nargs.size(), "fadd.size() != nargs.size()");
    
//    int ichosen = rand()%int(fadd.size());
    for (int i: RANDOM_ORDER1(int(fadd.size()))) {                    /* Filter to add */
//        if (!all && i != ichosen) { continue; }
        
        if (typeid(*fadd[i]) == typeid(FilterDownsamplePrefilter<real>)) {    /* Add downsample/upsample pair */
            int iup_min = 1, iup_max = int(f0->L.size());
            int idown_min = 0, idown_max = iup_max;
            if (!all) {
                iup_min = iup_max = 1+(rand()%f0->L.size());
                idown_min = idown_max = rand()%iup_min;
            }
            if (params.verbose >= 2) {
                printf("iup_min: %d, iup_max: %d\n", iup_min, iup_max);
                printf("idown_min: %d, idown_max: %d\n", idown_min, idown_max);
            }
            for (int iup: RANDOM_ORDER2(iup_min, iup_max+1)) {
                if (params.verbose >= 2) { printf("looping iup=%d, iup-1=%d, idown_max=%d, MIN(both) = %d\n", iup, iup-1, idown_max, MIN((iup-1), idown_max)); }
                if (iup < (int) f0->L.size() && f0->L[iup].inputs.size() > 1) { continue; }
                for (int idown: RANDOM_ORDER2(idown_min, MIN((iup-1), idown_max)+1)) {
                    if (params.verbose >= 2) { printf("looping iup=%d, idown=%d\n", iup, idown); }
                    ASSERT(idown < iup, "expected idown < iup");
                    
                    shared_ptr<FilterDAG<real> > f = f0->copy_sameclass();
                    if (iup < (int) f->L.size()) {
                        ASSERT(f->L[iup].inputs.size() == 1, "expected 1 input for L[iup]");
                        insert_dag(f, iup, shared_ptr<Filter<real> >(new FilterUpsamplePostfilter<real>()), f->L[iup].inputs);
                    } else {
                        vector<int> inputs_up({int(f->L.size())-1});
                        insert_dag(f, iup, shared_ptr<Filter<real> >(new FilterUpsamplePostfilter<real>()), inputs_up);
                    }
                    if (params.verbose >= 2) {
                        printf("After insert downsample:\n%s\n\n", f->str().c_str());
                    }
                    ASSERT(f->L[idown].inputs.size(), "expected L[idown] to have >= 1 inputs");
                    vector<int> inputs_down({f->L[idown].inputs[0]});
                    insert_dag(f, idown, shared_ptr<Filter<real> >(new FilterDownsamplePrefilter<real>()), inputs_down, true);
                    if (params.verbose >= 2) {
                        printf("Original:\n%s\n\n", f0->str().c_str());
                        printf("idown: %d, iup: %d\n", idown, iup);
                        printf("Result:\n%s\n\n", f->str().c_str());
                    }
                    ADD_NODE("add_mutate success for DownsamplePrefilter");
                }
            }
        } else {                                                      /* Add single node */
//            int jchosen = rand()%int(f0->L.size()+1);
            for (int jnode: RANDOM_ORDER1(int(f0->L.size())+1)) {      /* Node to add before */
//                if (!all && jnode != jchosen) { continue; }
                bool jlast = (jnode == f0->L.size());
                
                vector<int> jnode_inputs({int(f0->L.size())-1});
                if (!jlast) {
                    jnode_inputs = f0->L[jnode].inputs;
                    //printf("jnode_inputs from original node: %s\n", vector_to_str(jnode_inputs).c_str());
                } else {
                    //printf("jnode_inputs from size (%d): %s\n", int(f0->L.size()), vector_to_str(jnode_inputs).c_str());
                }
                ASSERT(jnode_inputs.size() >= 1, "jnode_inputs size is 0");
//                int kchosen = rand()%int(jnode_inputs.size());
//                for (int kinput: RANDOM_ORDER1(int(jnode_inputs.size()))) {
//                    if (!all && kinput != kchosen) { continue; }
                
                /* Input to new node */
                vector<vector<int> > inputs_lists;
                if (nargs[i] == 1) {
                    inputs_lists.push_back({jnode_inputs[0]});
                } else if (nargs[i] == 2) {
//                        int lchosen = (rand()%(jnode+1))-1;
                    for (int larg: RANDOM_ORDER2(-1, jnode)) {
//                            if (!all && larg != lchosen) { continue; }
                        inputs_lists.push_back({jnode_inputs[0], larg == -1 ? INPUT_SOURCE_IMAGE: larg});
                    }
                } else {
                    ASSERT(false, "3+ args not supported");
                }
                for (int lidx: RANDOM_ORDER1(int(inputs_lists.size()))) {
                    shared_ptr<FilterDAG<real> > f = f0->copy_sameclass();
                    //printf("Original f:\n%s\n\n", f->str().c_str());
                    //printf("insert_dag, jnode=%d, inputs_lists: %s\n", jnode, vector_to_str(inputs_lists[lidx]).c_str());
                    insert_dag(f, jnode, fadd[i]->copy(), inputs_lists[lidx]);
                    
                    ADD_NODE("add_mutate success");
                }
                    
//                }
            }
        }
    }
    return ans;
}

template<class real>
shared_ptr<FilterDAG<real> > get_subprogram_fragment(shared_ptr<FilterDAG<real> > f0) {
    if (params.verbose >= 1.5) {
        printf("get_subprogram_fragment begin, original program:\n%s\n\n", f0->str().c_str());
    }
    vector<double> input_size;
    
    if (f0->L.size() <= 1) { return f0->copy_sameclass(); }
    while (1) {
        int start = rand()%int(f0->L.size());
        if (params.verbose >= 2) {
            printf("get_subprogram iter, start=%d\n", start);
        }
        vector<int> dist = caller_distance(f0, start);
        int max_dist = *std::max_element(dist.begin()+start, dist.end());
        
        int chosen_dist = max_dist == 0 ? 0: (rand()%(max_dist+1));
        
        shared_ptr<FilterDAG<real> > ans = make_shared<FilterDAG<real> >();
        map<int, int> idx_map;
        for (int i = 0; i < (int) f0->L.size(); i++) {
            if (dist[i] >= 0 && dist[i] <= chosen_dist) {
                ASSERT(idx_map.count(i) == 0, "expected i not in idx_map");
                int prev_size = idx_map.size();
                idx_map[i] = prev_size;
                vector<int> inputs;
                for (int j = 0; j < (int) f0->L[i].inputs.size(); j++) {
                    int input = f0->L[i].inputs[j];
                    if (input == -1) { inputs.push_back(input); }
                    else {
                        if (!idx_map.count(input)) {
                            goto subprogram_failed;
                        } else {
                            inputs.push_back(idx_map[input]);
                        }
                    }
                }
                NodeDAG<real> node(f0->L[i].f->copy(), inputs);
                ans->L.push_back(node);
            }
        }
        
        if (ans->check_sizes(input_size)) {
            if (params.verbose >= 2) {
                printf("get_subprogram_fragment: %d\n", (int) ans->L.size());
            }
            return ans;
        }
        
        subprogram_failed:
        continue;
    }
}

template<class real>
shared_ptr<FilterDAG<real> > crossover_once(shared_ptr<FilterDAG<real> > a0, shared_ptr<FilterDAG<real> > b0) {
    if (params.verbose >= 1.5) {
        printf("crossover_once begin\n");
    }
    auto a = a0, b = b0;
    if (rand()%2 == 0) {
        std::swap(a, b);
    }
    if (rand()%2 == 0) {
        a = get_subprogram_fragment(a);
    } else {
        a = a->copy_sameclass();
    }
    b = get_subprogram_fragment(b);
    vector<int> valid_a({int(a->L.size())});              /* Indexes valid to insert 'b' graph before 'a' node */
    for (int i = 0; i < (int) a->L.size(); i++) {
        if (a->L[i].inputs.size() == 1) {
            valid_a.push_back(i);
        }
    }
    //printf("valid_a: %s\n", vector_to_str_int(valid_a).c_str());
    
    int idx_a = valid_a[rand()%int(valid_a.size())];
    int orig_input = int(a->L.size())-1;
    if (idx_a < (int) a->L.size()) {
        //printf("idx_a: %d, a->L[idx_a].inputs.size(): %d\n", idx_a, int(a->L[idx_a].inputs.size()));
        ASSERT(a->L[idx_a].inputs.size() == 1, "expected 1 input to a->L[idx_a]");
        orig_input = a->L[idx_a].inputs[0];
    }
    insert_subgraph(a, b, idx_a, orig_input);
    
    vector<double> input_size;
    a->check();
    //ASSERT(a->check_sizes(input_size), "a->check_sizes failed");
    if (!a->check_sizes(input_size)) { return shared_ptr<FilterDAG<real> >(); }
    if (params.verbose >= 1.5) {
        printf("crossover_once(%d, %d) => %d\n", int(a0->L.size()), int(b0->L.size()), int(a->L.size()));
    }
    return a;
}

template<class real>
vector<shared_ptr<FilterDAG<real> > > crossover(shared_ptr<FilterDAG<real> > a0, shared_ptr<FilterDAG<real> > b0, set<string> &seen) {
    if (params.verbose >= 1.5) {
        printf("crossover begin\n");
    }
    vector<shared_ptr<FilterDAG<real> > > ans;
    vector<double> input_size;
    bool all = false;
    
    for (int i = 0; i < params.ga_crossover_attempts; i++) {
        if (params.verbose >= 1.5) {
            printf("crossover iter %d\n", i);
        }
        auto f = crossover_once(a0, b0);
        if (f) {
            ADD_NODE("crossover success");
        }
    }
    if (params.verbose >= 1.5) {
        printf("crossover end no result\n");
    }
    return ans;
}

#endif
