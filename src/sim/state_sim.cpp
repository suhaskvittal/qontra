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
StateSimulator::reduce_record_by(uint64_t s) {
    record_offset -= s;
}

void
StateSimulator::shift_record_by(uint64_t offset) {
    for (uint64_t i = 0; i < record_offset; i++) {
        record_table[i].clear();
        record_table[i].swap_with(record_table[i + offset]);
    }
    record_offset -= offset;
}

void
StateSimulator::snapshot() {
    record_table_cpy = stim::simd_bit_table(record_table);
    record_offset_cpy = record_offset;
}

void
StateSimulator::rollback_at_trial(uint64_t t) {
    for (uint i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        record_table[i][t] = record_table_cpy[i][t];
    }
}

}   // qontra
