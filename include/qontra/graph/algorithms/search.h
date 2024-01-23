/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef GRAPH_ALGORITHMS_SEARCH_h
#define GRAPH_ALGORITHMS_SEARCH_h

#include "qontra/graph.h"

namespace qontra {
namespace graph {

// A function that can perform either BFS or DFS (code is essentially
// the same, except for the data structure used to store the vertices:
// queue vs stack, both of which can be modelled by a deque).
template <class V, class E, class FUNC>
void xfs(Graph<V, E>*, sptr<V> start, FUNC, bool dfs);

}   // graph
}   // qontra

#include "search.inl"

#endif  // GRAPH_ALGORITHMS_SEARCH_h
