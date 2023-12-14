/*
 *  author: Suhas Vittal
 *  date:   11 December 2023
 * */

#include "sim/memory_sim.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

namespace qontra {

using namespace graph;
using namespace lattice;

void
MemorySimulator::lrc_reset() {
    lrc_await_queue.clear();
    // Populate the await queue.
    for (uint s = 0; s < config.lrc_stride_size; s++) {
        for (uint pq : parity_qubits) {
            lrc_await_queue.push_back(pq);
        }
    }
    for (uint dq : data_qubits) lrc_await_queue.push_back(dq);
}

void
MemorySimulator::lrc_execute_lrcs_from_await_queue() {
    auto it = lrc_await_queue.begin();

    std::set<uint> qubits_in_use;

    std::deque<uint> lrc_await_queue_new_entries;
    // First go through and identify which parity qubits must be scheduled.
    std::vector<uint> data_qubits_awaiting_lrc;
    for (uint i = 0; i < k; i++) {
        uint q = lrc_await_queue.front();
        lrc_await_queue.pop_front();

        auto v = lattice_graph.get_vertex(q);
        if (v->qubit_type != vertex_t::type::data) {
            qubits_in_use.insert(q);

            lrc_await_queue_new_entries.push_back(q);
        } else {
            data_qubits_awaiting_lrc.push_back(q);
            lrc_await_queue_new_entries.push_back(q);
        }
    }
    // Now go through and deal with the data qubits.
    //
    // We need to solve a maximum matching problem here. Get
    // the available parity qubits.
    std::vector<uint> available_parity_qubits;
    for (uint pq : parity_qubits) {
        if (!qubits_in_use.count(pq)) available_parity_qubits.push_back(pq);
    }
    std::map<uint, uint> lrc_map =
        lrc_solve_maximum_matching(data_qubits_awaiting_lrc, available_parity_qubits);
    // Now, check if each data qubit in need of an LRC has been scheduled. If they could
    // not get an LRC, then remove them from lrc_await_queue_new_entries and place them
    // at the top of the lrc_await_queue.
    for (uint dq : data_qubits_awaiting_lrc) {
        if (!lrc_map.count(dq)) {
            for (auto it = lrc_await_queue_new_entries.begin();
                    it != lrc_await_queue_new_entries.end(); )
            {
                if (*it == dq)  it = lrc_await_queue_new_entries.erase(it);
                else            it++;
            }
            lrc_await_queue.push_front(dq);
        }
    }
    // Update lrc_await_queue.
    for (uint x : lrc_await_queue_new_entries) lrc_await_queue.push_back(x);
    // Finally do the LRC.
    lrc_measure_qubits(lrc_map);
}

void
MemorySimulator::lrc_optimal_oracle() {
    // Just need to identify which qubits are leaked.
}

void
MemorySimulator::lrc_measure_qubits(const std::map<uint, uint>& swap_map) {
    std::map<uint, uint> swap_map_inv;
    // Compute LRC operands and populate swap_map_inv.
    std::vector<uint> lrc_operands;
    for (auto& pair : swap_map) {
        lrc_operands.push_back(pair.first);
        lrc_operands.push_back(pair.second);

        swap_map_inv[pair.second] = pair.first;
    }
    std::vector<uint> lrc_operands_r(lrc_operands.rbegin(), lrc_operands.rend());
    // Perform any gates now.
    if (config.lrc_circuit == lrc_circuit_t::swap) {
        // Identify which qubits must be measured.
        std::vector<uint> measure_list;
        for (uint pq : parity_qubits) {
            if (swap_map_inv.count(pq)) {
                measure_list.push_back(swap_map_inv[pq]);
            } else {
                measure_list.push_back(pq);
            }
        }

        if (lrc_operands.size()) {
            elapsed_time += do_gate("cx", lrc_operands);
            inject_idling_error_negative(lrc_operands);
            elapsed_time += do_gate("cx", lrc_operands_r);
            inject_idling_error_negative(lrc_operands);
            elapsed_time += do_gate("cx", lrc_operands);
            inject_idling_error_negative(lrc_operands);
        }

        elapsed_time += do_measurement(measure_list);
        inject_idling_error_negative(measure_list);
        elapsed_time += do_gate("reset", measure_list);
        // We now have to update the measurement counters, as we measured some data qubits.
        for (auto pair : swap_map) {
            uint q1 = pair.first, q2 = pair.second;
            if (q1 == q2) continue;
            auto v1 = lattice_graph.get_vertex(q1);
            if (v1->qubit_type == vertex_t::type::data) {
                meas_ctr_map[q2] = meas_ctr_map[q1];
                meas_ctr_map.erase(q1);
            }
        }

        if (lrc_operands.size()) {
            elapsed_time += do_gate("cx", lrc_operands_r);
            inject_idling_error_negative(lrc_operands);
            elapsed_time += do_gate("cx", lrc_operands);
            inject_idling_error_negative(lrc_operands);
        }
    } else if (config.lrc_circuit == lrc_circuit_t::dqlr) {
        elapsed_time += do_measurement(parity_qubits);
        inject_idling_error_positive(data_qubits);
        elapsed_time += do_gate("reset", parity_qubits);

        elapsed_time += do_gate("liswap", lrc_operands_r);
        inject_idling_error_negative(lrc_operands);

        elapsed_time += do_gate("reset", parity_qubits);
    }
}

std::map<uint, uint>
MemorySimulator::lrc_solve_maximum_matching(
        const std::vector<uint>& avail_data_qubits,
        const std::vector<uint>& avail_parity_qubits)
{
    using namespace lemon;
    ListGraph gr;

    std::map<uint, uint> lrcs;

    // Create nodes and edges.
    std::map<uint, ListGraph::Node> qubit_to_graph_node;
    std::map<ListGraph::Node, uint> graph_node_to_qubit;

    for (uint pq : avail_parity_qubits) {
        ListGraph::Node x = gr.addNode();
        qubit_to_graph_node[pq] = x;
        graph_node_to_qubit[x] = pq;
    }

    for (uint dq : avail_data_qubits) {
        ListGraph::Node x = gr.addNode();
        qubit_to_graph_node[dq] = x;
        graph_node_to_qubit[x] = dq;

        auto dv = lattice_graph.get_vertex(dq);
        for (auto pv : lattice_graph.get_neighbors(dv)) {
            uint pq = pv->id;
            if (std::find(avail_parity_qubits.begin(),
                            avail_parity_qubits.end(),
                            pq) == avail_parity_qubits.end()) continue;
            ListGraph::Node y = qubit_to_graph_node[pq];
            gr.addEdge(x, y);
        }
    }
    // Run maximum matching algorithm.
    MaxMatching m(gr);
    m.run();
    for (uint dq : avail_data_qubits) {
        auto data_node = qubit_to_graph_node[dq];
        auto parity_node = m.mate(data_node);
        if (parity_node == INVALID) continue;
        uint pq = graph_node_to_qubit[parity_node];
        lrcs[dq] = pq;
    }
    return lrcs;
}

}   // qontra
