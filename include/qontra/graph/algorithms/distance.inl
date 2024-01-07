/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

#include <limits>
#include <queue>
#include <vector>

#include <math.h>

namespace qontra {
namespace graph {

// Performs Dijkstra's algorithm given an edge weight function (ewf_t).
template <class V, class E, class W_FUNC> void
dijkstra(
    Graph<V, E>* graph, 
    sptr<V> src,
    std::map<sptr<V>, fp_t>& distances,
    std::map<sptr<V>, sptr<V>>& predecessors,
    W_FUNC edge_w_func)
{
    typedef struct { sptr<V> v; fp_t s; } pqv_t;
    struct cmp {
        bool operator()(const pqv_t& v1, const pqv_t& v2) {
            return v1.s > v2.s;
        }
    };

    std::map<sptr<V>, pqv_t> v2pv;
    std::priority_queue<pqv_t, std::vector<pqv_t>, cmp> queue;
    for (sptr<V> v : graph->get_vertices()) {
        if (v == src)   distances[v] = 0;
        else            distances[v] = std::numeric_limits<fp_t>::max();
        predecessors[v] = v;

        pqv_t pv = {v, distances[v]};
        queue.push(pv);
        v2pv[v] = pv;
    }

    while (!queue.empty()) {
        pqv_t pv = queue.top();
        sptr<V> v = pv.v;
        queue.pop();
        if (fabsl(pv.s - distances[v]) > 1e-8) continue;   // This entry is outdated.

        std::vector<sptr<V>> adj = graph->get_neighbors(v);
        for (sptr<V> w : adj) {
            sptr<E> e = graph->get_edge(v, w);
            fp_t new_dist = distances[v] + edge_w_func(e);
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
template <class V, class E, class DATA, class W_FUNC, class DATA_FUNC> DistanceMatrix<V, DATA>
create_distance_matrix(Graph<V, E>* graph, W_FUNC edge_w_func, DATA_FUNC cb) {
    DistanceMatrix<V, DATA> mat;
    std::vector<sptr<V>> vertices = graph->get_vertices();
    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<V> src = vertices[i];
        std::map<sptr<V>, sptr<V>> pred;
        std::map<sptr<V>, fp_t> dist;
        dijkstra(graph, src, dist, pred, edge_w_func);
        for (size_t j = 0; j < vertices.size(); j++) {
            if (i == j) continue;
            sptr<V> dst = vertices[j];
            DATA x = cb(src, dst, dist, pred);
            vtils::tlm_put(mat, src, dst, x);
        }
    }
    return mat;
}

}   // graph
}   // qontra
