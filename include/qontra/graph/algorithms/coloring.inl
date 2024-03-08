/*
 *  author: Suhas Vittal
 *  date:   25 February 2024
 * */

#include <algorithm>
#include <deque>

namespace qontra {
namespace graph {

template <class V, class E> int
k_coloring_greedy(Graph<V, E>* graph, std::map<sptr<V>, int>& color_map, size_t bfs_seed) {
    std::vector<sptr<V>> vertices = graph->get_vertices();
    if (vertices.empty()) return -1;
    // Solve by BFS.
    std::set<sptr<V>> visited;
    int max_color = -1;
    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<V> st = vertices.at( (i + bfs_seed) % vertices.size() );
        std::deque<sptr<V>> bfs{ st };
        while (bfs.size()) {
            sptr<V> v = bfs.front();
            bfs.pop_front();
            if (visited.count(v)) continue;
            // Try and color v by examining its neighbors.
            std::set<int> used;
            for (sptr<V> w : graph->get_neighbors(v)) {
                bfs.push_back(w);
                if (color_map.count(w)) {
                    used.insert(color_map.at(w));
                }
            }
            int c = 0;
            while (true) {
                if (used.count(c)) {
                    c++;
                } else {
                    color_map[v] = c;
                    break;
                }
            }
            max_color = std::max(max_color, c);
            visited.insert(v);
        }
        if (visited.size() == vertices.size()) break;
    }
    return max_color;
}

template <class V, class E> int
k_coloring_rlf(Graph<V, E>* graph, std::map<sptr<V>, int>& color_map) {
    std::vector<sptr<V>> vertices = graph->get_vertices();

    // Track vertex deletions.
    std::map<sptr<V>, size_t> degree_offset_map;

    int c = -1;
    while (vertices.size()) {
        c++;
        std::set<sptr<V>> color_set;
        // Find largest vertex amongst vertices.
        auto max_it = std::max_element(vertices.begin(), vertices.end(),
                [&] (sptr<V> x, sptr<V> y)
                {
                    return graph->get_degree(x) - degree_offset_map[x] 
                            < graph->get_degree(y) - degree_offset_map[y];
                });
        color_set.insert(*max_it);

        std::set<sptr<V>> adjacent_to_color_set;
        for (sptr<V> w : graph->get_neighbors(*max_it)) {
            if (!color_map.count(w)) {
                adjacent_to_color_set.insert(w);
                degree_offset_map[w]++;
            }
        }
        vertices.erase(max_it);
        while (true) {
            // Now, find vertices not adjacent to any vertices in color_set that have a maximal
            // amount of neighbors adjacent to vertices in color_set.
            int neighbors_adj_to_color_set = -1;
            typename std::vector<sptr<V>>::iterator best_it;
            for (auto it = vertices.begin(); it != vertices.end(); it++) {
                if (adjacent_to_color_set.count(*it)) continue;
                int k = 0;
                for (sptr<V> w : graph->get_neighbors(*it)) {
                    if (!color_map.count(w)) {
                        k += adjacent_to_color_set.count(w);
                    }
                }
                if (k > neighbors_adj_to_color_set
                    || (k == neighbors_adj_to_color_set
                        && graph->get_degree(*it) - degree_offset_map[*it] < 
                            graph->get_degree(*best_it) - degree_offset_map[*best_it])
                    )
                {
                    neighbors_adj_to_color_set = k;
                    best_it = it;
                }
            }
            if (neighbors_adj_to_color_set >= 0) {
                color_set.insert(*best_it);
                for (sptr<V> w : graph->get_neighbors(*best_it)) {
                    if (!color_map.count(w)) {
                        adjacent_to_color_set.insert(w);
                        degree_offset_map[w]++;
                    }
                }
                vertices.erase(best_it);
            } else {
                break;
            }
        }
        // Color all vertices in color_set as c.
        for (sptr<V> v : color_set) color_map[v] = c;
    }
    return c;
}

}   // graph
}   // qontra
