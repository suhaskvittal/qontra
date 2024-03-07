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
    W_FUNC edge_w_func,
    sptr<V> target)
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
    /*
    typedef std::pair<int, fp_t> pqv_t;
    struct cmp {
        bool operator()(const pqv_t& v1, const pqv_t& v2) const {
            return v1.second > v2.second;
        }
    };

    std::priority_queue<pqv_t, std::vector<pqv_t>, cmp> queue;

    const size_t n = graph->n();
    std::vector<fp_t> _dist(n, std::numeric_limits<fp_t>::max());
    std::vector<int> _pred(n, -1);

//  const auto& enum_map = graph->get_enumeration_map();
    auto vertices = graph->get_vertices();
    std::map<sptr<V>, int> enum_map;
    for (size_t i =0; i < vertices.size(); i++) enum_map[vertices[i]] = i;

    _dist[enum_map.at(src)] = 0;
    queue.emplace(enum_map.at(src), 0);

    while (!queue.empty()) {
        auto [ i, d ] = queue.top();
        queue.pop();
        if (d > _dist[i]) continue;   // This entry is outdated.
        
        sptr<V> v = vertices.at(i);
        if (target != nullptr && v == target) break;

        for (sptr<V> w : graph->get_neighbors(v)) {
            int j = enum_map.at(w);
            sptr<E> e = graph->get_edge(v, w);
            fp_t new_dist = d + edge_w_func(e);
            if (new_dist < _dist[j]) {
                _dist[j] = new_dist;
                _pred[j] = i;
                queue.emplace(j, new_dist);
            }
        }
    }
    // Update distances, predecesors.
    for (size_t i = 0; i < n; i++) {
        sptr<V> v = vertices.at(i);
        distances[v] = _dist[i];
        
        int j = _pred[i];
        predecessors[v] = j < 0 ? v : vertices.at(j);
    }
    */
}

template <class V, class E, class W_FUNC> void
floyd_warshall(
    Graph<V, E>* graph, 
    vtils::TwoLevelMap<sptr<V>, sptr<V>, fp_t>& distances,
    vtils::TwoLevelMap<sptr<V>, sptr<V>, sptr<V>>& predecessors,
    W_FUNC edge_w_func)
{
    const size_t n = graph->n();
    // Represent as vectors first as this is more efficient. Translate to the
    // two-level map later.
    std::vector<std::vector<fp_t>> _dist(n, std::vector<fp_t>(n, std::numeric_limits<fp_t>::max()));
    std::vector<std::vector<int>> _pred(n);

    auto vertices = graph->get_vertices();
    for (size_t i = 0; i < n; i++) {
        _dist[i][i] = 0;
        _pred[i] = std::vector<int>(n, i);
    }
    const auto& enum_map = graph->get_enumeration_map();
    for (sptr<E> e : graph->get_edges()) {
        sptr<V> v = std::reinterpret_pointer_cast<V>(e->src),
                w = std::reinterpret_pointer_cast<V>(e->dst);
        size_t i = enum_map.at(v),
               j = enum_map.at(w);
        _dist[i][j] = edge_w_func(e);
        _dist[j][i] = _dist[i][j];
    }
    // Compute distances.
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            fp_t d_ik = _dist[i][k];
            for (size_t j = i+1; j < n; j++) {
                fp_t& d_ij = _dist[i][j],
                    & d_ji = _dist[j][i];
                fp_t d_kj = _dist[k][j];
                if (d_ik != std::numeric_limits<fp_t>::max() 
                    && d_kj != std::numeric_limits<fp_t>::max()
                    && d_ij > d_ik + d_kj) 
                {
                    d_ij = d_ik + d_kj;
                    d_ji = d_ij;
                    _pred[i][j] = _pred[k][j];
                    _pred[j][i] = _pred[k][i];
                }
            }
        }
    }
    // Move to two-level-map.
    for (size_t i = 0; i < n; i++) {
        sptr<V> v = vertices.at(i);
        distances[v][v] = 0;
        predecessors[v][v] = v;
        for (size_t j = i+1; j < n; j++) {
            sptr<V> w = vertices.at(j);

            distances[v][w] = _dist[i][j];
            distances[w][v] = _dist[j][i];

            sptr<V> u1 = vertices.at(_pred[i][j]),
                    u2 = vertices.at(_pred[j][i]);
            predecessors[v][w] = u1;
            predecessors[w][v] = u2;
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
