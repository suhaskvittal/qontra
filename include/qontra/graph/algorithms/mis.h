/* author: Suhas Vittal
 *  date:   14 June 2023
 * */

#ifndef GRAPH_ALGORITHMS_MIS_h
#define GRAPH_ALGORITHMS_MIS_h

#include "qontra/graph.h"

#include <set>

namespace qontra {
namespace graph {

// Algorithms for computing the maximum independent set of a graph:
//
// mis_approx_greedy achieves a (1 / (D+1)) approximation ratio, where D is the
// maximum degree of the graph. Furthermore, the minimum size of the output
// is (n / (2(d+1))) where d is the mean degree of the graph and n is the
// number of nodes.

template <class V, class E>
std::set<sptr<V>> mis_approx_greedy(Graph<V, E>* graph);

}   // graph
}   // qontra

#include "mis.inl"

#endif
