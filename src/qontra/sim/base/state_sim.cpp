/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "qontra/sim/base/state_sim.h"

namespace qontra {

size_t G_RECORD_SPACE_SIZE = 8192;

StateSimulator::StateSimulator(uint64_t n, uint64_t _max_shots)
    :record_table(G_RECORD_SPACE_SIZE, _max_shots),
    lock_table(n, _max_shots),
    shots(_max_shots),
    record_table_cpy(G_RECORD_SPACE_SIZE, _max_shots),
    lock_table_cpy(n, _max_shots),
    rng(0),
    n_qubits(n),
    max_shots(_max_shots)
{
    reset_sim();
}

StateSimulator::StateSimulator(const StateSimulator& other)
    :record_table(other.record_table),
    lock_table(other.lock_table),
    shots(other.shots),
    record_table_cpy(other.record_table_cpy),
    lock_table_cpy(other.lock_table_cpy),
    rng(other.rng),
    n_qubits(other.n_qubits),
    max_shots(other.max_shots)
{}

}   // qontra
