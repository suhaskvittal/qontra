/*
 *  author: Suhas Vittal
 *  date:   20 June 2023
 * */

#include "sim/control_sim.h"

namespace qontra {

using namespace experiments;

void
ControlSimulator::IF() {
    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if_stall[t] = decoder_busy[t].not_zero();
        if (!if_stall[t])   if_pc[t].u64[0] = pc[t].u64[0];
        if_id_valid[t] = 1;
    }
}

void
ControlSimulator::ID() {
    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        auto inst = program[if_pc[t].u64[0]];
    }
}

}   // qontra
