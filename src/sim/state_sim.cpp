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
        if (i < record_offset - offset) {
            record_table[i].swap_with(record_table[i + offset]);
        } else {
            record_table[i].clear();
        }
    }
    record_offset -= offset;
}

}   // qontra
