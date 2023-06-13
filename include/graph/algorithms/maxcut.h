/*
 *  author: Suhas Vittal
 *  date:   13 June 2023
 * */

#ifndef GRAPH_ALGORITHMS_MAXCUT_h
#define GRAPH_ALGORITHMS_MAXCUT_h

#include "graph/graph.h"

#include <set>
#include <vector>

namespace qontra {
namespace graph {
namespace cut {

template <class V_t, int k>
using Partition = std::array<std::set<V_t>, k>;

namespace maxcut {

template <class V_t, class E_t> Partition<V_t, 2>
maxcut_approx_greedy(Graph<V_t, E_t>* graph, ewf_t edge_w_func) {
    auto vertices = graph->get_vertices();

    Partition<V_t, 2>  partition;
    partition.fill(std::set<V_t>());

    for (auto v : vertices) {
        fp_t bias = 0;   // If bias < 0, then put v in the right set.
        for (auto w : graph->get_neighbors(v)) {
            fp_t x = edge_w_func(v, w);
            if (partition[0].count(w))          bias -= x;
            else if (partition[1].count(w))     bias += x;
        }
        partition[bias < 0.0].insert(v);
    }
    return partition;
}

}

}   // cut
}   // graph
}   // qontra

#endif  // GRAPH_ALGORITHMS_MAXCUT_h
