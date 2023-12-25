/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "sim/state_sim.h"

namespace qontra {

namespace statesim {
    uint64_t G_RECORD_SPACE_SIZE = 4096;
}

void
StateSimulator::shift_record_by(uint64_t offset) {
    for (uint64_t i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        if (i < offset) record_table[i].clear();
        else            record_table[i].swap_with(record_table[i-offset]);
    }
}

void
StateSimulator::snapshot() {
    record_table_cpy = stim::simd_bit_table(record_table);
    lock_table_cpy = stim::simd_bit_table(lock_table);
}

void
StateSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    for (uint i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        copy_where(record_table_cpy[i], record_table[i], pred);
    }
    for (uint i = 0; i < n_qubits; i++) {
        copy_where(lock_table_cpy[i], lock_table[i], pred);
    }
}

}   // qontra
