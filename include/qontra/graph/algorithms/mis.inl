/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

#include <algorithm>

namespace qontra {
namespace graph {

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
