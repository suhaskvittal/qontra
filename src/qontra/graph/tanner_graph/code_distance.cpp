/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "qontra/graph/tanner_graph.h"

#include <stim.h>

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
        push_back_operator(circuit, 
                            gr->get_neighbors(op),
                            op->qubit_type == vertex_t::type::xparity);
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

}   // graph
}   // qontra
