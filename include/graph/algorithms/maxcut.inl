/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

namespace qontra {
namespace graph {

template <class V, class E, class W_FUNC> partition_t<V>
maxcut_approx_greedy(Graph<V, E>* graph, W_FUNC edge_w_func) {
    std::vector<sptr<V>> vertices = graph->get_vertices();

    partition_t<V>  partition;
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
