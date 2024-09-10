/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef GRAPH_ALGORITHMS_DISTANCE_h
#define GRAPH_ALGORITHMS_DISTANCE_h

#include "qontra/graph.h"

#include <vtils/two_level_map.h>

namespace qontra {
namespace graph {

template <class V, class DATA>
using DistanceMatrix = vtils::TwoLevelMap<sptr<V>, sptr<V>, DATA>;

// Performs Dijkstra's algorithm given an edge weight function. If target is null,
// then the algorithm computes the distance to all vertices from src.
template <class V, class E, class W_FUNC>
void dijkstra(
        Graph<V, E>*,
        sptr<V> src,
        std::unordered_map<sptr<V>, fp_t>& dist,
        std::unordered_map<sptr<V>, sptr<V>>& pred,
        W_FUNC,
        sptr<V> target=nullptr);

template <class V, class E, class W_FUNC>
void floyd_warshall(
        Graph<V, E>*,
        vtils::TwoLevelMap<sptr<V>, sptr<V>, fp_t>& dist,
        vtils::TwoLevelMap<sptr<V>, sptr<V>, sptr<V>>& prev,
        W_FUNC);

template <class V, class E, class DATA, class W_FUNC, class DATA_FUNC>
DistanceMatrix<V, DATA> create_distance_matrix(Graph<V, E>*, W_FUNC, DATA_FUNC);

}   // graph
}   // qontra

#include "distance.inl"

#endif  // GRAPH_ALGORITHMS_DISTANCE_h
