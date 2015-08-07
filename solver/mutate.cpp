#include "mutate.h"

void upsample_downsample_dfs(int start, vector<bool> &visited, map<int, set<int> > &edges, bool &terminate_early) {
    ASSERT((unsigned) start < visited.size(), "start out of bounds");
    if (!visited[start]) {
        visited[start] = true;
        if (visited[visited.size()-1]) {
            terminate_early = true;
            return;
        }
        if (edges.count(start)) {
            for (int neighbor: edges[start]) {
                upsample_downsample_dfs(neighbor, visited, edges, terminate_early);
                if (terminate_early) { return; }
            }
        }
    }
}

vector<int> random_order(int lo, int hi, bool is_random) {
    int n = hi - lo;
    if (n <= 0) { return {}; }
    vector<int> L(n);
    for (int i = 0; i < n; i++) {
        L[i] = i + lo;
    }
    if (is_random) {
        std::random_shuffle(L.begin(), L.end());
    }
    return L;
}

vector<int> random_order(int n, bool is_random) {
    return random_order(0, n, is_random);
}
