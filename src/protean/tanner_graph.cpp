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
    if (v->qubit_type == vertex_t::DATA)    return pred;

    std::set<vertex_t*> visited;
    auto v_adj = get_neighbors(v);
    for (auto d : v_adj) {
        auto d_adj = get_neighbors(d);
        for (auto w : d_adj) {
            auto w_adj = get_neighbors(w);
            if (v_adj.size() != w_adj.size() && !visited.count(w) && is_subset_of(w_adj, v_adj)) {
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

    std::sort(v1_adj.begin(), v1_adj.end());
    std::sort(v2_adj.begin(), v2_adj.end());

    bool either_precedes = is_subset_of(v1_adj, v2_adj) || is_subset_of(v2_adj, v1_adj); 
    if (either_precedes) return nullptr;

    std::vector<vertex_t*> induced_adj;
    std::set_intersection(v1_adj.begin(), v1_adj.end(), v2_adj.begin(), v2_adj.end(),
                            std::back_inserter(induced_adj));
    if (induced_adj.size() < 2)    return nullptr;
    // If there is any intersection, create a vertex and add it and the edges to the graph.
    auto w = new vertex_t;
    w->qubit_type = vertex_t::GAUGE;
    w->id = (induced_gauge_index++) | INDUCED_GAUGE_INDEX_FLAG | (w->qubit_type << 30);
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

namespace graph {
namespace io {

void
update_tanner_graph(protean::TannerGraph& graph, std::string line) {
    using namespace protean;

    size_t ssi;   // Variable for calculating substring indices.

    if (line.size() == 0)   return; // Nothing to be done.
    bool is_x_check = line[0] == 'X';
    // Get the check id.
    ssi = line.find(",", 1);
    if (ssi == std::string::npos)   return;
    uint check_id = std::stoi(line.substr(1, ssi));

    // Create check vertex.
    tanner::vertex_t* check_v = new tanner::vertex_t;
    check_v->qubit_type = is_x_check ? tanner::vertex_t::XPARITY : tanner::vertex_t::ZPARITY;
    check_v->id = check_id | (check_v->qubit_type << 30);
    if (!graph.add_vertex(check_v)) {
        delete check_v;
        return;
    }
    // Now, parse the data qubits involved in the check.
    uint pssi = ssi+1;
    while (true) {
        if ((ssi = line.find(",", pssi)) == std::string::npos) {
            ssi = line.size();
        }
        uint dqid = std::stoi(line.substr(pssi, ssi));
        auto data_v = graph.get_vertex(dqid);
        if (data_v == nullptr) {
            data_v = new tanner::vertex_t;
            data_v->id = dqid;
            data_v->qubit_type = tanner::vertex_t::DATA;
            graph.add_vertex(data_v);
        }
        // Add edge.
        tanner::edge_t* e = new tanner::edge_t;
        e->src = check_v;
        e->dst = data_v;
        graph.add_edge(e);  // This should never fail.
        pssi = ssi+1;
        if (ssi == line.size()) break;
    } 
}

}   // io
}   // graph

}   // qontra
