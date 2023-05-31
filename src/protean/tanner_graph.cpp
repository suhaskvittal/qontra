/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "protean/tanner_graph.h"

namespace qontra {
namespace protean {

std::vector<tanner::vertex_t*>
TannerGraph::get_predecessors(tanner::vertex_t* v) {
    std::vector<tanner::vertex_t*> pred;
    std::set<tanned::vertex_t*> visited;

    auto v_adj = adjacency_list(v);
    for (auto d : v_adj) {
        auto d_adj = adjacency_list(d);
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

tanner::vertex_t*
TannerGraph::induce_predecessor(tanner::vertex_t* v1, tanner::vertex_t* v2) {
    auto v1_adj = adjacency_list(v1);
    auto v2_adj = adjacency_list(v2);

    std::vector<tanner::vertex_t*> induced_adj;
    std::set_intersection(v1_adj.begin(), v1_adj.end(), v2_adj.begin(), v2_adj.end(),
                            std::back_inserter(induced_adj));
    if (induced_adj.size() == 0)    return nullptr;
    // If there is any intersection, create a vertex and add it and the edges to the graph.
    auto w = new tanner::vertex_t;
    w->id = (induced_gauge_index++) | INDUCED_GAUGE_INDEX_FLAG;
    add_vertex(w);
    for (auto u : induced_adj) {
        auto e =  new tanner::edge_t;
        e->src = v;
        e->dst = u;
        add_edge(e);
    }

    return w;
}

}   // protean
}   // qontra
