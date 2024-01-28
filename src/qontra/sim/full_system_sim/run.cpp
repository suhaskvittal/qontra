/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "qontra/sim/full_system_sim.h"

namespace qontra {

FullSystemSimulator::FullSystemSimulator()
    :register_file(1, 1),
    elapsed_time(0),
    shot_time_delta_map(),
    syndrome_table(1, 1),
    observable_table(1, 1),
    meas_ctr(0),
    meas_offset(0),
    event_offset(0),
    max_event_written_to(0),
    max_obs_written_to(0),
    sample_circuit(),
    is_recording_sim_instructions(false),
    n_qubits(1),
    current_shots(0),
    base_sim(nullptr),
    subroutine_map(),
    register_file_cpy(1, 1),
    syndrome_table_cpy(1, 1),
    observable_table_cpy(1, 1)
{}

void
FullSystemSimulator::run_batch(const qes::Program<>& program, uint64_t shots) {
    base_sim->shots = shots;
    current_shots = shots;
    // Reset all structures.
    register_file.clear();
    
    elapsed_time = 0.0;
    shot_time_delta_map.clear();

    syndrome_table.clear();
    observable_table.clear();
    
    meas_ctr = 0;
    meas_offset = 0;
    event_offset = 0;
    max_event_written_to = 0;
    max_obs_written_to = 0;

    execute_routine(program);
}

}   // qontra
