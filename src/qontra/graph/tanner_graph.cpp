/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "qontra/graph/tanner_graph.h"
#include "qontra/graph/algorithms/coloring.h"

#include <stim.h>

#include <limits>

namespace qontra {
namespace graph {

using namespace tanner;

void push_back_operator(stim::Circuit& circuit, std::vector<sptr<vertex_t>> op, bool is_x) {
    std::vector<uint32_t> targets;
    for (sptr<vertex_t> dq : op) {
        uint32_t q = static_cast<uint32_t>(dq->id);
        if (is_x) {
            targets.push_back(q | stim::TARGET_PAULI_X_BIT);
        } else {
            targets.push_back(q | stim::TARGET_PAULI_Z_BIT);
        }
        targets.push_back(stim::TARGET_COMBINER);
    }
    targets.pop_back();
    circuit.safe_append_u("MPP", targets);
}

void measure_stabilizers(stim::Circuit& circuit, const TannerGraph* gr) {
    for (sptr<vertex_t> op : gr->get_checks()) {
        push_back_operator(circuit, gr->get_neighbors(op), op->qubit_type == vertex_t::type::xparity);
    }
}

void measure_logical_qubits(stim::Circuit& circuit, const TannerGraph* gr, bool is_x) {
    size_t obsno = 0;
    auto obs_list = gr->get_obs(is_x);
    for (auto& obs : obs_list) {
        push_back_operator(circuit, obs, is_x);
        std::vector<uint32_t> obs_targets{ 1 | stim::TARGET_RECORD_BIT };
        circuit.safe_append_ua("OBSERVABLE_INCLUDE", obs_targets, static_cast<double>(obsno));
        obsno++;
    }
}

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
    for (size_t s = 0; s < checks.size() && (best_coloring > 2 || best_coloring < 0); s++) {
        std::map<sptr<vertex_t>, int> cm;
        int c = k_coloring_greedy(gr.get(), cm, s);
        if (c < best_coloring) {
            check_color_map = std::move(cm);
            cm.clear();
            best_coloring = c;
        }
    }
    return best_coloring;
}

int
TannerGraph::compute_code_distance(bool is_x) const {
    stim::Circuit mem;
    measure_logical_qubits(mem, this, is_x);
    measure_stabilizers(mem, this);
    // Inject errors.
    std::vector<uint32_t> targets;
    for (sptr<vertex_t> v : get_vertices_by_type(vertex_t::type::data)) {
        targets.push_back(static_cast<uint32_t>(v->id));
    }
    mem.safe_append_ua("DEPOLARIZE1", targets, 1e-3);
    measure_stabilizers(mem, this);
    // Add detection events.
    const size_t n_checks = get_checks().size();
    for (size_t i = 0; i < n_checks; i++) {
        targets = std::vector<uint32_t>{ 
                        static_cast<uint32_t>(i+1) | stim::TARGET_RECORD_BIT,
                        static_cast<uint32_t>(i+1+n_checks) | stim::TARGET_RECORD_BIT
                    };
        mem.safe_append_ua("DETECTOR", targets, static_cast<double>(i));
    }
    measure_logical_qubits(mem, this, is_x);
    // Analyze errors and get distance.
    auto dem =
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
                mem,
                false,  // decompose errors
                true,   // fold loops
                false,  // allow gauge detectors (non-deterministic detectors)
                0.0,    // approx disjoint errors threshold
                false,  // ignore decomposition failures
                false
            );
    auto errors = stim::find_undetectable_logical_error(dem, 100, 10, false);
    return errors.count_errors();
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
