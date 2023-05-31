/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "graph/coupling_graph.h"

namespace qontra {
namespace graph {

bool
CouplingGraph::test_planarity_after_add(const std::vector<coupling::edge_t*>& edges) {
    if (!is_planar) return false;   // Adding edges can only make something nonplanar. No need
                                    // to check if the graph is already nonplanar.
    bool curr_planarity = is_planar;

    // Temporarily update the graph
    for (auto e : edges)    Graph::add_edge(e);
    check_planarity();
    bool will_be_planar = is_planar;
    for (auto e : edges)    Graph::delete_edge(e);

    is_planar = curr_planarity;
    return will_be_planar;
}

bool
CouplingGraph::test_planarity_after_delete(const std::vector<coupling::vertex_t*>& vertices) {
    if (is_planar)  return true;    // Deleting vertices can only make something planar. No need
                                    // to check if the graph is already planar.
    bool curr_planarity = is_planar;

    std::set<coupling::edge_t*> removed_edges;
    for (auto v : vertices) {
        auto adj = adjacency_list(v);   // As a vertex deletion also removes edges, we need
                                        // to figure out which edges will be removed.
        for (auto w : adj) {
            auto e = Graph::get_edge(v, w);
            removed_edges.insert(e);
        }
        Graph::delete_vertex(v);
    }
    check_planarity();
    bool will_be_planar = is_planar;
    for (auto v : vertices)         Graph::add_vertex(v);
    for (auto e : removed_edges)    Graph::add_edge(e);

    is_planar = curr_planarity;
    return will_be_planar;
}

bool
CouplingGraph::test_planarity_after_delete(const std::vector<coupling::edge_t*>& edges) {
    if (is_planar)  return true;    // See rationale in the above function.
    bool curr_planarity = is_planar;
    
    for (auto e : edges)    Graph::delete_edge(e);
    check_planarity();
    bool will_be_planar = is_planar;
    for (auto e : edges)    Graph::add_edge(e);

    is_planar = curr_planarity;
    return will_be_planar;
}

void
CouplingGraph::check_planarity() {
    if (!track_planarity) return;
    if (edges.size() > 3*vertices.size() - 5) {
        is_planar = false;
        return;
    }
    // We'll use LEMON for this.
    lemon::ListGraph planar_graph;
    std::map<coupling::vertex_t*, lemon::ListGraph::Node> lemon_node_map;
    for (auto v : vertices) {
        auto x = planar_graph.add_node();
        lemon_node_map[v] = x;
    }
    for (auto e : edges) {
        auto x = lemon_node_map[e->src];
        auto y = lemon_node_map[e->dst];
        planar_graph.add_arc(x, y);
    }
    is_planar = lemon::checkPlanarity(planar_graph);
}

}   // graph
}   // qontra
