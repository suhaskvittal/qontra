/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

#include "qontra/experiments.h"

namespace qontra {

inline size_t
get_register_index(std::string r) {
    const size_t REGISTER_SPECIAL_START = 32;
    // Check if the register is special, otherwise it is straightforward to get the index.
    if (r == "rbrk") {
        return REGISTER_SPECIAL_START + 0;
    } else {
        int id = std::stoi(r.substr(1));
        return static_cast<size_t>(id);
    }
}

template <class SIM> void
FullSystemSimulator::run_program(const qes::Program<>& program, uint64_t shots) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    // Execute program simulation in batches:
    const uint64_t local_shots = shots / world_size + static_cast<int>(world_rank==0)*(shots % world_size);
    base_sim->set_seed(G_BASE_SEED + world_rank);

    // Set up all structures.
    is_recording_stim_instructions = true;
    n_qubits = get_number_of_qubits(program);
    base_sim = uptr(new SIM(n_qubits, G_SHOTS_PER_BATCH));

    register_file = stim::simd_bit_table<SIMD_WIDTH>(n_registers, G_SHOTS_PER_BATCH);
    syndrome_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);
    observable_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);

    uint64_t shots_remaining = local_shots;
    while (shots_remaining) {
        const uint64_t shots_this_batch = shots_remaining < local_shots ? shots_remaining : G_SHOTS_PER_BATCH;
        run_batch(program, shots_this_batch);
        // Save stats.
        is_recording_stim_instructions = false;
        shots_remaining -= shots_this_batch;
    }
}

inline void
FullSystemSimulator::load_subroutine(std::string name, const qes::Program<>& program) {
    subroutine_map[name] = program;
}

inline void
FullSystemSimulator::execute_routine(const qes::Program<>& program) {
    program_status_t status;
    while (status.pc < program.size()) {
        const qes::Instruction<>& instruction = program.at(status.pc);
        // Check if the instruction is requesting the execution of the
        // microcode.
        if (instruction.get_name() == "call") {
            // Switch into the subroutine.
            std::string subroutine_name = instruction.get<std::string>(1);
            execute_routine(subroutine_map.at(subroutine_name));
            status.pc++;
        } else {
            read_next_instruction(current_program, main_status);
            if (status.return_if_waiting_trials.popcnt() == current_shots) break;
        }
    }
}

inline void
FullSystemSimulator::snapshot() {
    base_sim->snapshot();

    register_file_cpy = register_file;
    syndrome_table_cpy = syndrome_table;
    observable_table_cpy = observable_table;
}

inline void
FullSystemSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    base_sim->rollback_where(pred);
    copy_where(register_file_cpy, register_file, pred);
    copy_where(syndrome_table_cpy, syndrome_table, pred);
    copy_where(observable_table_cpy, observable_table, pred);
}

}   // qontra
