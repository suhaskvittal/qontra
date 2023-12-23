/*
 *  author: Suhas Vittal
 *  date:   11 December 2023
 * */

#include "sim/memory_sim.h"

namespace qontra {

using namespace graph;
using namespace lattice;

void
MemorySimulator::eraser_initialize() {
    for (uint i = 0; i < parity_qubits.size(); i++) {
        eraser_qubit_to_syndrome_buffer_index_map[parity_qubits[i]] = i;
    }
    // Setup swap lookup table (depth = 2).
    for (uint dq : data_qubits) {
        uint earliest_cx_time = std::numeric_limits<uint>::max();
        uint latest_cx_time = 0;
        sptr<vertex_t> pv_earliest = nullptr;
        sptr<vertex_t> pv_latest = nullptr;

        auto dv = lattice_graph.get_vertex(dq);
        for (auto pv : lattice_graph.get_neighbors(dv)) {
            auto e = lattice_graph.get_edge(dv, pv);
            uint t = e->cx_time;
            if (t < earliest_cx_time) {
                earliest_cx_time = t;
                pv_earliest = pv;
            }
            if (t > latest_cx_time) {
                latest_cx_time = t;
                pv_latest = pv;
            }
        }
        eraser_swap_lookup_table[dq] = std::array<uint, 2>();
        eraser_swap_lookup_table[dq][0] = pv_earliest->id;
        eraser_swap_lookup_table[dq][1] = pv_latest->id;
    }
    eraser_recently_scheduled_qubits_table =
        stim::simd_bit_table(experiments::G_SHOTS_PER_BATCH, n_qubits);
    eraser_syndrome_buffer = 
        stim::simd_bit_table(parity_qubits.size(), experiments::G_SHOTS_PER_BATCH);
    eraser_syndrome_buffer.clear();
}

void
MemorySimulator::eraser_reset() {
    eraser_recently_scheduled_qubits_table.clear();
    eraser_syndrome_buffer.clear();
}

stim::simd_bits
MemorySimulator::eraser_execute_lrcs() {
    // Examine the syndromes from the most recent round.
    //
    // (1) If any data qubit is adjacent to #neighbors/2 flipped parity checks, then
    //      schedule it for an LRC.
    // (2) Do not use any recently scheduled qubits.
    stim::simd_bit_table eraser_syndrome_buffer_tr = eraser_syndrome_buffer.transposed();

    stim::simd_bits shots_with_leakage(sim->shots);
    shots_with_leakage.clear();
    for (uint64_t t = 0; t < sim->shots; t++) {
        if (eraser_syndrome_buffer_tr[t].popcnt() == 0) continue;

        std::map<uint, uint> lrc_map;
        std::set<uint> qubits_in_use;
        for (uint dq : data_qubits) {
            if (eraser_recently_scheduled_qubits_table[t][dq]) continue;
            auto dv = lattice_graph.get_vertex(dq);
            // Count number of flipped neighboring checks.
            const uint dg = lattice_graph.get_degree(dv);
            uint flipped_checks = 0;
            for (auto pv : lattice_graph.get_neighbors(dv)) {
                uint pq = pv->id;
                uint k = eraser_qubit_to_syndrome_buffer_index_map[pq];
                flipped_checks += eraser_syndrome_buffer_tr[t][k];
            }
            if (flipped_checks < dg/2) continue;
            // Try and schedule data qubit for an LRC.
            for (uint pq : eraser_swap_lookup_table[dq]) {
                bool in_use = qubits_in_use.count(pq)
                                || eraser_recently_scheduled_qubits_table[t][pq];
                if (!in_use) {
                    lrc_map[dq] = pq;
                    qubits_in_use.insert(dq);
                    qubits_in_use.insert(pq);
                    break;
                }
            }
        }
        eraser_recently_scheduled_qubits_table[t].clear();
        // If we have not scheduled any LRCs, move on.
        if (lrc_map.empty()) continue;
        // Otherwise, update the recently_scheduled_qubits_table and perform LRCs.
        for (uint q : qubits_in_use) eraser_recently_scheduled_qubits_table[t][q] = 1;
        lrc_measure_qubits(lrc_map, (int64_t)t);
        shots_with_leakage[t] = 1;
    }
    return shots_with_leakage;
}

void
MemorySimulator::eraser_update_syndrome_buffer(const std::map<uint, uint64_t>& prev_meas_ctr_map) {
    eraser_syndrome_buffer.clear();
    for (uint pq : parity_qubits) {
        uint k = eraser_qubit_to_syndrome_buffer_index_map[pq];
        if (prev_meas_ctr_map.count(pq)) {
            uint64_t pmt = prev_meas_ctr_map.at(pq);
            eraser_syndrome_buffer[k] ^= sim->record_table[pmt];
        }
        uint64_t mt = meas_ctr_map[pq];
        eraser_syndrome_buffer[k] ^= sim->record_table[mt];
    }
}

}   // qontra
