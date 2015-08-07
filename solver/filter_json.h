
#ifndef _filter_json_h
#define _filter_json_h

#include <json/json.h>
#include "util.h"

#undef ASSERT
#define ASSERT(x, msg) ASSERT2(x, msg)

int name_to_index(const string &a, int n);

class ProgramNameSort { public:
    int n;
    ProgramNameSort(int n_) :n(n_) { }
    bool operator() (const string &a, const string &b) {
        int ai = name_to_index(a, n);
        int bi = name_to_index(b, n);
        return ai < bi;
    }
};

template<class real>
Array<real> parse_array2d(Json::Value a) {
    ASSERT(a.isArray(), "parse_array2d expected array argument");
    ASSERT(a[0u].isArray(), "parse_array2d expected array of array argument");
    int h = a.size(), w = a[0u].size();
    vector<int> sizes({h, w});
    Array<real> ans(sizes);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            ans(y, x) = a[y][x].asDouble();
        }
    }
    return ans;
}

template<class real>
Array<real> parse_array1d(Json::Value a) {
    ASSERT(a.isArray(), "parse_array1d expected array argument");
    int n = a.size();
    vector<int> sizes({n});
    Array<real> ans(sizes);
    for (int i = 0; i < n; i++) {
        ans(i) = a[i].asDouble();
    }
    return ans;
}

template<class real>
shared_ptr<FilterSparse<real> > read_program(Json::Value program, Json::Value mask, bool has_mask) {
    auto ans = make_shared<FilterDAG<real> >();
    ASSERT(program.isObject(), "expected program (DAG) to be object");
    auto keys = program.getMemberNames();
    int n = keys.size()-1;
    sort(keys.begin(), keys.end(), ProgramNameSort(n));
    ASSERT(keys.size(), "expected program (DAG) to be non-empty");
    for (int i = 0; i < (int) keys.size(); i++) {
        Json::Value sub(program[keys[i]]);
        ASSERT(sub.isArray(), "expected DAG element to be array");
        ASSERT(sub.size() >= 1, "expected DAG element to have size >= 1");
        string filter_type = sub[0u].asString();
        if (filter_type == "ImageParam") {
        } else if (filter_type == "FIRFilter") {
            ASSERT(sub.size() == 3, "expected FIRFilter to have size 3");
            int input = name_to_index(sub[1u].asString(), n);
            auto A = parse_array2d<real>(sub[2u]);
            shared_ptr<Filter<real> > f;
            if (A.height() == 1) {
                f = shared_ptr<Filter<real> >(new FilterFIRH<real>(A));
            } else if (A.width() == 1) {
                f = shared_ptr<Filter<real> >(new FilterFIRV<real>(A));
            } else {
                f = shared_ptr<Filter<real> >(new FilterFIR<real>(A));
            }
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else if (filter_type == "IIRFilter") {
            ASSERT(sub.size() == 6, "expected IIRFilter to have size 6");
            int input = name_to_index(sub[1u].asString(), n);
            Array<real> K = parse_array1d<real>(sub[2u]);
            Array<real> F = parse_array1d<real>(sub[3u]);
            int dx = sub[4u].asInt();
            int dy = sub[5u].asInt();
            shared_ptr<Filter<real> > f;
            if      (dx ==  1 && dy ==  0) { f = shared_ptr<Filter<real> >(new FilterIIR<real>(K, F)); }
            else if (dx == -1 && dy ==  0) { f = shared_ptr<Filter<real> >(new FilterIIRMinusX<real>(K, F)); }
            else if (dx ==  0 && dy == -1) { f = shared_ptr<Filter<real> >(new FilterIIRMinusY<real>(K, F)); }
            else if (dx ==  0 && dy ==  1) { f = shared_ptr<Filter<real> >(new FilterIIRPlusY<real>(K, F)); }
            else {
                printf("Warning: got IIRFilter in exotic direction %d, %d\n", dx, dy);
                f = shared_ptr<Filter<real> >(new FilterIIR<real>(K, F, dx, dy));
            }
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else if (filter_type == "FilterAdd") {
            ASSERT(sub.size() >= 3, "expected FilterAdd to have size >= 3");
            vector<int> inputL;
            for (unsigned int idx = 1u; idx < sub.size(); idx++) {
                inputL.push_back(name_to_index(sub[idx].asString(), n));
            }
            ASSERT2(inputL.size() >= 2, "expected >= 2 args for FilterAdd");
            auto f = make_shared<FilterAdd<real> >();
            ans->L.push_back(NodeDAG<real>(f, inputL));
        } else if (filter_type == "FilterUpsample") {
            ASSERT(sub.size() == 2, "expected FilterUpsample to have size 2");
            int input = name_to_index(sub[1u].asString(), n);
            auto f = make_shared<FilterUpsample<real> >();
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else if (filter_type == "FilterDownsample") {
            ASSERT(sub.size() == 2, "expected FilterDownsample to have size 2");
            int input = name_to_index(sub[1u].asString(), n);
            auto f = make_shared<FilterDownsample<real> >();
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else if (filter_type == "FilterUpsamplePostfilter") {
            ASSERT(sub.size() == 2, "expected FilterUpsamplePostfilter to have size 2");
            int input = name_to_index(sub[1u].asString(), n);
            auto f = make_shared<FilterUpsamplePostfilter<real> >();
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else if (filter_type == "FilterDownsamplePrefilter") {
            ASSERT(sub.size() == 2, "expected FilterDownsamplePrefilter to have size 2");
            int input = name_to_index(sub[1u].asString(), n);
            auto f = make_shared<FilterDownsamplePrefilter<real> >();
            ans->L.push_back(NodeDAG<real>(f, {input}));
        } else {
            ASSERT(0, "read_program: parsed unknown filter_type");
        }
    }
    
    auto ans_sparse = make_shared<FilterSparse<real> >(ans);
    if (has_mask) {
        ASSERT(mask.isArray(), "expected mask to be array");
        if (ans_sparse->mask.size() != mask.size()) {
            fprintf(stderr, "ans_sparse mask size: %d, mask size: %d\n", (int) ans_sparse->mask.size(), (int) mask.size());
            ASSERT(false, "expected ans_sparse mask size == mask size");
        }
        for (unsigned i = 0; i < mask.size(); i++) {
            ans_sparse->mask[i] = int_to_masktype(mask[i].asInt());
        }
    } else {
        ans_sparse->mask_from_coeffs();
    }
    return ans_sparse;
}

template<class real>
vector<shared_ptr<FilterSparse<real> > > json_to_filters(string filename, vector<double> *scaleL=NULL) {
    if (scaleL) { scaleL->clear(); }
    vector<shared_ptr<FilterSparse<real> > > ans;
    
    string file_str = read_binary_file(filename);
    Json::Value root;
    Json::Reader reader;
    bool parsingSuccessful = reader.parse(file_str, root);
    if (!parsingSuccessful) {
        fprintf(stderr, "could not parse json file %s\n", filename.c_str());
        exit(1);
    }
    
    if (!root.isArray()) {
        fprintf(stderr, "Error: json root is not array. Loading %s\n", filename.c_str());
    }
    ASSERT(root.isArray(), "json root is not array");
    for (int i = 0; i < (int) root.size(); i++) {
        Json::Value sub = root[i];
        ASSERT(sub.isObject(), "json array element is not object (dict)");
        ASSERT(sub.isMember("name"), "sub does not contain field 'name'");
        string name = sub["name"].asString();
        ASSERT(sub.isMember("program"), "sub does not contain field 'program'");
        Json::Value program = sub["program"];
        bool has_mask = false;
        Json::Value mask;
        if (sub.isMember("mask")) {
//        ASSERT(sub.isMember("mask"), "sub does not contain field 'mask'");
            mask = sub["mask"];
            has_mask = true;
        }
        if (scaleL && sub.isMember("scale")) {
            scaleL->push_back(sub["scale"].asDouble());
        }
        ans.push_back(read_program<real>(program, mask, has_mask));
    }
    return ans;
}

template<class real>
vector<shared_ptr<FilterDAG<real> > > json_to_filters_dag(string filename) {
    vector<shared_ptr<FilterDAG<real> > > ans;
    vector<shared_ptr<FilterSparse<real> > > L = json_to_filters<real>(filename);
    ASSERT2(L.size(), "expected to load 1 or more filters from json");
    for (int i = 0; i < (int) L.size(); i++) {
        ASSERT2(typeid(*L[i]->filter) == typeid(FilterDAG<real>), "json did not load FilterSparse(FilterDAG())");
        shared_ptr<FilterDAG<real> > dag = dynamic_pointer_cast<FilterDAG<real> >(L[i]->filter);
        ans.push_back(dag);
    }
    return ans;
}

#endif

