/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

#include <deque>
#include <set>

namespace qontra {
namespace graph {

// A function that can perform either BFS or DFS (code is essentially
// the same, except for the data structure used to store the vertices:
// queue vs stack, both of which can be modelled by a deque).
template <class V, class E, class FUNC> void
xfs(Graph<V, E>* graph, sptr<V> start, FUNC cb, bool dfs) {
    std::deque<sptr<V>> dq;
    std::set<sptr<V>> visited;
    dq.push_back(start);

    while (dq.size()) {
        sptr<V> v;
        if (dfs) {
            v = dq.back();
            dq.pop_back();
        } else {
            v = dq.front();
            dq.pop_front();
        }

        for (auto w : graph->get_neighbors(v)) {
            if (visited.count(w))   continue;
            cb(v, w);
            dq.push_back(w);
        }
        visited.insert(v);
    }
}

}   // graph
}   // qontra
