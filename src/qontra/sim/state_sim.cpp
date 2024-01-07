/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "qontra/sim/state_sim.h"

namespace qontra {

namespace statesim {
    size_t G_RECORD_SPACE_SIZE = 8192;
}

StateSimulator::StateSimulator(uint n, uint64_t max_shots)
    :n_qubits(n),
    max_shots(max_shots),
    record_table(statesim::G_RECORD_SPACE_SIZE, max_shots),
    lock_table(n, max_shots),
    record_table_cpy(statesim::G_RECORD_SPACE_SIZE, max_shots),
    lock_table_cpy(n, max_shots),
    rng(0)
{
    reset_sim();
}

StateSimulator::StateSimulator(const StateSimulator& other)
    :n_qubits(other.n_qubits),
    max_shots(other.max_shots),
    record_table(other.record_table),
    lock_table(other.lock_table),
    record_table_cpy(other.record_table_cpy),
    lock_table_cpy(other.lock_table_cpy),
    rng(other.rng)
{}

}   // qontra
