/*
 *  author: Suhas Vittal
 *  date:   13 June 2023
 * */
#ifndef GRAPH_ALGORITHMS_MAXCUT_h
#define GRAPH_ALGORITHMS_MAXCUT_h

#include "qontra/graph.h"

#include <array>
#include <set>

namespace qontra {
namespace graph {

template <class V, size_t K=2>
using partition_t = std::array<std::set<V*>, K>;

// Algorithms for computing the max-cut of a graph:
//
// approx_greedy achieves a 0.5 approximation ratio.

template <class V, class E, class W_FUNC>
partition_t<V> maxcut_approx_greedy(Graph<V, E>*, W_FUNC);

}   // graph
}   // qontra

#include "maxcut.inl"

#endif  // GRAPH_ALGORITHMS_MAXCUT_h
