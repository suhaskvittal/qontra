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
}

void
StateSimulator::rollback_at_trial(uint64_t t) {
    for (uint i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        record_table[i][t] = record_table_cpy[i][t];
    }
}

}   // qontra
