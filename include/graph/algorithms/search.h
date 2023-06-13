/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef GRAPH_ALGORITHMS_SEARCH_h
#define GRAPH_ALGORITHMS_SEARCH_h

#include "graph/graph.h"

namespace qontra {
namespace graph {
namespace search {  // For BFS/DFS

template <class V_t>
using callback_t = std::function<void(V_t*, V_t*)>; // This callback is
                                                    // called on a 
                                                    // successful traversal
                                                    // from a vertex
                                                    // to its unvisited
                                                    // neighbor


// A function that can perform either BFS or DFS (code is essentially
// the same, except for the data structure used to store the vertices:
// queue vs stack, both of which can be modelled by a deque).
template <class V_t, class E_t> void
xfs(Graph<V_t, E_t>* graph, V_t* start, callback_t<V_t> cb, bool dfs) {
    std::deque<V_t*> dq;
    std::set<V_t*> visited;
    dq.push_back(start);

    while (dq.size()) {
        V_t* v;
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

}   // search
}   // graph
}   // qontra

#endif  // GRAPH_ALGORITHMS_SEARCH_h
