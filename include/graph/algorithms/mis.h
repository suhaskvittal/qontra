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
namespace mis {

// Algorithms for computing the maximum independent set of a graph:
//
// approx_greedy achieves a (1 / (D+1)) approximation ratio, where D is the
// maximum degree of the graph. Furthermore, the minimum size of the output
// is (n / (2(d+1))) where d is the mean degree of the graph and n is the
// number of nodes.

template <class V_t, class E_t> std::set<V_t*>
approx_greedy(Graph<V_t, E_t>* graph) {
    auto vertices = graph->get_vertices();
    std::set<V_t*> marked;

    auto cmp = [&] (V_t* x, V_t* y) 
    {
        return graph->get_degree(x) < graph->get_degree(y);   
    };
    std::sort(vertices.begin(), vertices.end(), cmp);

    std::set<V_t*> indep;
    for (auto v : vertices) {
        if (marked.count(v))    continue;
        // Otherwise, add v to the IS and mark its neighbors.
        indep.insert(v);
        for (auto w : graph->get_neighbors(v))  marked.insert(w);
    }
    return indep;
}

}   // mis
}   // graph
}   // qontra

#endif
