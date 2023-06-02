/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "protean/tanner_graph.h"

namespace qontra {
namespace protean {

using namespace tanner;

std::vector<vertex_t*>
TannerGraph::get_predecessors(vertex_t* v) {
    std::vector<vertex_t*> pred;
    std::set<vertex_t*> visited;

    auto v_adj = get_neighbors(v);
    for (auto d : v_adj) {
        auto d_adj = get_neighbors(d);
        for (auto w : d_adj) {
            if (v_adj.size() != d_adj.size() && !visited.count(w)
                && std::includes(d_adj.begin(), d_adj.end(), v_adj.begin(), v_adj.end()))
            {
                pred.push_back(w);
                visited.insert(w);
            }
        }
    }
    return pred;
}

vertex_t*
TannerGraph::induce_predecessor(vertex_t* v1, vertex_t* v2) {
    auto v1_adj = get_neighbors(v1);
    auto v2_adj = get_neighbors(v2);

    bool either_precedes = is_subset_of(v1_adj, v2_adj) || is_subset_of(v2_adj, v1_adj); 
    if (either_precedes) return nullptr;

    std::vector<vertex_t*> induced_adj;
    std::set_intersection(v1_adj.begin(), v1_adj.end(), v2_adj.begin(), v2_adj.end(),
                            std::back_inserter(induced_adj));
    if (induced_adj.size() == 0)    return nullptr;
    // If there is any intersection, create a vertex and add it and the edges to the graph.
    auto w = new vertex_t;
    w->id = (induced_gauge_index++) | INDUCED_GAUGE_INDEX_FLAG;
    add_vertex(w);
    for (auto u : induced_adj) {
        auto e = new edge_t;
        e->src = (void*)w;
        e->dst = (void*)u;
        add_edge(e);
    }

    return w;
}

}   // protean
}   // qontra
