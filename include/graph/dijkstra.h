/*
 *  author: Suhas Vittal
 *  date:   18 December 2022
 * */

#ifndef GRAPH_DIJKSTRA_h
#define GRAPH_DIJKSTRA_h

#include "defs.h"

#include <functional>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <vector>

namespace qrc {
namespace graph {
    
template <class V> 
struct DijkstraResult {
    std::vector<V*> path;
    fp_t distance;
};

template <class V>
using PathTable = std::map<std::pair<V*, V*>, DijkstraResult<V>>;

template <class G, class V> void
dijkstra(G& graph,
        V * src, 
        std::map<V*, fp_t>& distances,
        std::map<V*, V*>& predecessors,
        std::function<fp_t(G&, V*, V*)>& edge_weight_func)
{
    typedef std::pair<V*, fp_t> PQVertex;

    struct DijkstraCmp {
        bool operator()(const PQVertex& v1, const PQVertex& v2) {
            return v1.second > v2.second; 
        };
    };
    std::map<V*, PQVertex> v2pv;
    std::priority_queue<PQVertex, std::vector<PQVertex>, DijkstraCmp> queue;
    for (V * v : graph.vertices()) {
        if (v == src) {
            distances[v] = 0; 
        } else {
            distances[v] = std::numeric_limits<fp_t>::max();
        } 
        predecessors[v] = v;

        PQVertex pv = std::make_pair(v, distances[v]);
        queue.push(pv);
        v2pv[v] = pv;
    }

    std::set<V*> visited;
    while (!queue.empty()) {
        PQVertex pvi = queue.top(); 
        auto vi = pvi.first;
        queue.pop();
        if (pvi.second != distances[vi]) {
            continue; 
        }

        auto adj_list = graph.adjacency_list(vi);
        for (V * vj : adj_list) {
            if (visited.count(vj)) {
                continue; 
            }
            fp_t new_dist = distances[vi] + edge_weight_func(graph, vi, vj);
            if (new_dist < distances[vj]) {
                distances[vj] = new_dist; 
                predecessors[vj] = vi;
                // Insert new entry into the priority queue.
                PQVertex pvj = std::make_pair(vj, new_dist);
                queue.push(pvj);
            }
        }
        visited.insert(vi);
    } 
}

template <class V> void
update_path_table(
        PathTable<V>& path_table, 
        V * src,
        V * dst,
        std::map<V*, fp_t>& distances,
        std::map<V*, V*>& predecessors)
{
    // Compute path.
    V * curr = dst;
    fp_t distance = distances[dst];
    std::vector<V*> path;

    while (curr != src) {
        if (curr == predecessors[curr]) {
            distance = std::numeric_limits<fp_t>::max();
            path.clear();
            goto failed_to_find_path;
        }
        path.push_back(curr);
        curr = predecessors[curr];
    }
    path.push_back(curr);
failed_to_find_path:
    // Build result.
    DijkstraResult<V> res = {path, distance};
    path_table[std::make_pair(src, dst)] = res;
    path_table[std::make_pair(dst, src)] = res;
}

}   // graph
}   // qrc

#endif
