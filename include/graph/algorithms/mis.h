/*
 *  author: Suhas Vittal
 *  date:   14 June 2023
 * */

#ifndef GRAPH_ALGORITHMS_MIS_h
#define GRAPH_ALGORITHMS_MIS_h

#include "graph/graph.h"

#include <algorithm>
#include <set>

namespace qontra {
namespace graph {

// Algorithms for computing the maximum independent set of a graph:
//
// mis_approx_greedy achieves a (1 / (D+1)) approximation ratio, where D is the
// maximum degree of the graph. Furthermore, the minimum size of the output
// is (n / (2(d+1))) where d is the mean degree of the graph and n is the
// number of nodes.

template <class V, class E> std::set<sptr<V>>
mis_approx_greedy(Graph<V, E>* graph) {
    std::vector<sptr<V>> vertices = graph->get_vertices();
    std::set<sptr<V>> marked;

    auto cmp = [&] (sptr<V> x, sptr<V> y) 
    {
        return graph->get_degree(x) < graph->get_degree(y);   
    };
    std::sort(vertices.begin(), vertices.end(), cmp);

    std::set<sptr<V>> indep;
    for (sptr<V> v : vertices) {
        if (marked.count(v))    continue;
        // Otherwise, add v to the IS and mark its neighbors.
        indep.insert(v);
        for (auto w : graph->get_neighbors(v))  marked.insert(w);
    }
    return indep;
}

}   // graph
}   // qontra

#endif
