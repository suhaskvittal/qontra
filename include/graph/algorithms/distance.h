/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef GRAPH_ALGORITHMS_DISTANCE_h
#define GRAPH_ALGORITHMS_DISTANCE_h

#include "graph/graph.h"

#include <math.h>

namespace qontra {
namespace graph {
namespace distance {

template <class V_t, class data_t>
using DistanceMatrix = TwoLevelMap<V_t*, V_t*, data_t>;

// The callback below is for creating a distance matrix over multiple
// calls of dijkstra's. The function should return some data for a
// matrix entry given (1, 2) two vertices indexing into the matrix,
// (3) the distances between the first input and every other vertex,
// (4) the predecessor entry for each vertex in a path from the 
// first input to that vertex.
template <class V_t, class data_t>
using callback_t = std::function<data_t(V_t*, V_t*, 
                        const std::map<V_t*, fp_t>&, const std::map<V_t*, V_t*>&)>;

// Performs Dijkstra's algorithm given an edge weight function (ewf_t).
template <class V_t, class E_t> void
dijkstra(
    Graph<V_t, E_t>* graph, 
    V_t* src,
    std::map<V_t*, fp_t>& distances,
    std::map<V_t*, V_t*>& predecessors,
    ewf_t<V_t> edge_w_func)
{
    typedef struct { V_t* v; fp_t s; } pqv_t;
    struct cmp {
        bool operator()(const pqv_t& v1, const pqv_t& v2) {
            return v1.s > v2.s;
        }
    };

    std::map<V_t*, pqv_t> v2pv;
    std::priority_queue<pqv_t, std::vector<pqv_t>, cmp> queue;
    for (V_t* v : graph->get_vertices()) {
        if (v == src)   distances[v] = 0;
        else            distances[v] = std::numeric_limits<fp_t>::max();
        predecessors[v] = v;

        pqv_t pv = {v, distances[v]};
        queue.push(pv);
        v2pv[v] = pv;
    }

    while (!queue.empty()) {
        pqv_t pv = queue.top();
        auto v = pv.v;
        queue.pop();
        if (fabsl(pv.s - distances[v]) > 1e-8) continue;   // This entry is outdated.

        auto adj = graph->get_neighbors(v);
        for (V_t* w : adj) {
            fp_t new_dist = distances[v] + edge_w_func(v, w);
            if (new_dist < distances[w]) {
                distances[w] = new_dist;
                predecessors[w] = v;
                pqv_t pw = {w, new_dist};
                queue.push(pw);
            }
        }
    }
}

// Creates a distance matrix by calling dijkstra's multiple times.
// Each entry is populated according to the callback (defined above in the
// header).
template <class V_t, class E_t, class data_t> DistanceMatrix<V_t, data_t>
create_distance_matrix(
    Graph<V_t, E_t>* graph, 
    ewf_t<V_t> edge_w_func,
    callback_t<V_t, data_t> cb) 
{
    DistanceMatrix<V_t, data_t> mat;
    auto vertices = graph->get_vertices();
    for (uint i = 0; i < vertices.size(); i++) {
        V_t* src = vertices[i];
        std::map<V_t*, V_t*> pred;
        std::map<V_t*, fp_t> dist;
        dijkstra(graph, src, dist, pred, edge_w_func);
        for (uint j = 0; j < vertices.size(); j++) {
            if (i == j) continue;
            V_t* dst = vertices[j];
            data_t x = cb(src, dst, dist, pred);
            tlm::put(mat, src, dst, x);
        }
    }
    return mat;
}

}   // distance
}   // graph
}   // qontra

#endif  // GRAPH_ALGORITHMS_DISTANCE_h
