/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "qontra/graph/tanner_graph.h"

#include "qontra/graph/algorithms/coloring.h"

#include <limits>

namespace qontra {
namespace graph {

using namespace tanner;

int
TannerGraph::update_check_color_map(std::map<sptr<vertex_t>, int>& check_color_map, bool use_x_checks) const {
    const std::vector<sptr<vertex_t>>& checks = use_x_checks ? xparity_checks : zparity_checks;
    // Create an interaction graph for these checks.
    uptr<Graph<vertex_t, base::edge_t>> gr = std::make_unique<Graph<vertex_t, base::edge_t>>();
    for (sptr<vertex_t> v : checks) gr->add_vertex(v);
    // Add edges if two vertices share neighbors in the Tanner graph.
    for (size_t i = 0; i < checks.size(); i++) {
        sptr<vertex_t> v = checks.at(i);
        for (size_t j = i+1; j < checks.size(); j++) {
            sptr<vertex_t> w = checks.at(j);
            if (get_common_neighbors({v, w}).size()) {
                gr->make_and_add_edge(v, w);
            }
        }
    }
    int best_coloring = k_coloring_rlf(gr.get(), check_color_map);
    for (size_t s = 0; s < checks.size() && best_coloring > 2; s++) {
        std::map<sptr<vertex_t>, int> cm;
        int c = k_coloring_greedy(gr.get(), cm, s);
        if (c < best_coloring) {
            check_color_map = std::move(cm);
            best_coloring = c;
        }
    }
    return best_coloring;
}

namespace io {

void
update_tanner_graph(TannerGraph& graph, std::string line) {
    size_t ssi;   // Variable for calculating substring indices.

    if (line.size() == 0)   return; // Nothing to be done.
    bool is_obs = line[0] == 'O';
    bool is_x_obs = line[1] == 'X';
    bool is_x_check = line[0] == 'X';
    // Get the check id.
    int check_id;
    ssi = line.find(",", 1);
    if (ssi == std::string::npos)   return;
    if (is_obs) {
        check_id = -1;
    } else {
        check_id = std::stoi(line.substr(1, ssi));
    }

    // Create check vertex.
    sptr<tanner::vertex_t> check_v;
    if (check_id >= 0) {
        const uint64_t flag = is_x_check ? VERTEX_ID_XPARITY_FLAG : VERTEX_ID_ZPARITY_FLAG;

        check_v = std::make_shared<tanner::vertex_t>();
        check_v->qubit_type = is_x_check ? tanner::vertex_t::type::xparity : tanner::vertex_t::type::zparity;
        check_v->id = static_cast<uint64_t>(check_id) | flag;
        if (!graph.add_vertex(check_v)) {
            return;
        }
    } else {
        check_v = nullptr;
    }
    // Now, parse the data qubits involved in the check.
    TannerGraph::obs_t obs_qubits;
    size_t pssi = ssi+1;
    while (true) {
        if ((ssi = line.find(",", pssi)) == std::string::npos) {
            ssi = line.size();
        }
        uint64_t dqid = static_cast<uint64_t>(std::stoi(line.substr(pssi, ssi-pssi))) | VERTEX_ID_DATA_FLAG;
        auto data_v = graph.get_vertex(dqid);
        if (data_v == nullptr) {
            data_v = std::make_shared<tanner::vertex_t>();
            data_v->id = dqid;
            data_v->qubit_type = tanner::vertex_t::type::data;
            graph.add_vertex(data_v);
        }
        // If this is a check, add an edge between the check and data qubit.
        // Otherwise, add the data qubit to the observable list.
        if (is_obs) {
            obs_qubits.push_back(data_v);
        } else {
            auto e = std::make_shared<tanner::edge_t>();
            e->src = std::static_pointer_cast<void>(check_v);
            e->dst = std::static_pointer_cast<void>(data_v);
            e->is_undirected = true;
            graph.add_edge(e); // This should never fail.
        }
        pssi = ssi+1;
        if (ssi == line.size()) break;
    }

    if (is_obs) {
        if (is_x_obs) {
            graph.x_obs_list.push_back(obs_qubits);
        } else {
            graph.z_obs_list.push_back(obs_qubits);
        }
    }
}

}   // io

}   // graph
}   // qontra
