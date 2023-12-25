/*
 *  author: Suhas Vittal
 *  date:   13 June 2023
 * */
#ifndef GRAPH_ALGORITHMS_MAXCUT_h
#define GRAPH_ALGORITHMS_MAXCUT_h

#include "graph/graph.h"

#include <array>
#include <set>

namespace qontra {
namespace graph {

template <class V, size_t K>
using partition_t = std::array<std::set<V*>, K>;

// Algorithms for computing the max-cut of a graph:
//
// approx_greedy achieves a 0.5 approximation ratio.

template <class V, class E, class W_FUNC> partition_t<V, 2>
maxcut_approx_greedy(Graph<V, E>* graph, W_FUNC edge_w_func) {
    std::vector<sptr<V>> vertices = graph->get_vertices();

    partition_t<V, 2>  partition;
    partition.fill(std::set<V>());

    for (sptr<V> v : vertices) {
        fp_t bias = 0;   // If bias < 0, then put v in the right set.
        for (sptr<V> w : graph->get_neighbors(v)) {
            sptr<E> e = graph->get_edge(v, w);
            fp_t x = edge_w_func(e);
            if (partition[0].count(w))          bias -= x;
            else if (partition[1].count(w))     bias += x;
        }
        partition[bias < 0].insert(v);
    }
    return partition;
}

}

}   // graph
}   // qontra

#endif  // GRAPH_ALGORITHMS_MAXCUT_h
