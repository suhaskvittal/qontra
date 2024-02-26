/*
 *  author: Suhas Vittal
 *  date:   25 February 2024
 * */

#ifndef GRAPH_ALGORITHMS_COLORING_h
#define GRAPH_ALGORITHMS_COLORING_h

#include "qontra/graph.h"

#include <map>

namespace qontra {
namespace graph {

template <class V, class E>
int k_coloring_greedy(Graph<V, E>*, std::map<sptr<V>, int>&);

template <class V, class E>
int k_coloring_rlf(Graph<V, E>*, std::map<sptr<V>, int>&);

}   // graph
}   // qontra

#include "coloring.inl"

#endif  // GRAPH_ALGORITHMS_COLORING_h
