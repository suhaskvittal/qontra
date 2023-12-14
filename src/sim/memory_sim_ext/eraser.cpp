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
    // Setup swap lookup table (depth = 2).
    for (auto dv : data_qubits) {
        uint earliest_cx_time = std::numeric_limits<uint>::max();
        uint latest_cx_time = 0;
        sptr<vertex_t> pv_earliest = nullptr;
        sptr<vertex_t> pv_latest = nullptr;
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
        eraser_swap_lookup_table[i] = std::array<uint, 2>();
        eraser_swap_lookup_table[i][0] = pv_earliest->id;
        eraser_swap_lookup_table[i][1] = pv_latest->id;
    }
    eraser_syndrome_buffer = stim::simd_bit_table(parity_qubits.size(), experiments::G_SHOTS_PER_BATCH);
    eraser_syndrome_buffer.clear();
}

void
MemorySimulator::eraser_reset() {
    eraser_recently_scheduled_qubits.clear();
}

std::map<uint, uint>
MemorySimulator::eraser_make_lrc_decisions() {
    std::map<uint, uint> lrc_schedule;
    // Examine the syndromes from the most recent round.
    //
    // (1) If any data qubit is adjacent to #neighbors/2 parity checks, then
    //      schedule it for an LRC.
    // (2) Do not use any recently scheduled qubits.
    for (auto dv : data_qubits) {
    }
    return lrc_schedule;
}

}   // qontra
