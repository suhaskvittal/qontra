/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

#include "qontra/experiments.h"

#include <vtils/filesystem.h>

namespace qontra {

inline size_t
get_register_index(std::string r) {
    const size_t REGISTER_SPECIAL_START = 32;
    // Check if the register is special, otherwise it is straightforward to get the index.
    if (r == "$rbrk") {
        return REGISTER_SPECIAL_START + 0;
    } else if (r == "$reoz") {
        return REGISTER_SPECIAL_START + 1;
    } else {
        int id = std::stoi(r.substr(2));
        return static_cast<size_t>(id);
    }
}

template <class SIM> histogram_t<uint64_t>
FullSystemSimulator::run_program(const qes::Program<>& program, uint64_t shots) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    // Execute program simulation in batches:
    const uint64_t local_shots = shots / world_size + static_cast<int>(world_rank==0)*(shots % world_size);

    // Set up all structures.
    is_recording_stim_instructions = true;
    n_qubits = get_number_of_qubits(program);
    base_sim = uptr<SIM>(new SIM(n_qubits, G_SHOTS_PER_BATCH));

    base_sim->set_seed(G_BASE_SEED + world_rank);

    register_file = stim::simd_bit_table<SIMD_WIDTH>(config.n_registers, G_SHOTS_PER_BATCH);
    syndrome_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);
    observable_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);

    shot_histogram.clear();

    // Create the stats files if they do not exist:
    using namespace vtils;
    if (!file_exists(config.syndrome_output_folder)) safe_create_directory(config.syndrome_output_folder);
    if (!file_exists(config.stim_output_file)) {
        safe_create_directory(get_parent_directory(config.stim_output_file.c_str()));
    }
    if (!file_exists(config.data_output_file)) {
        safe_create_directory(get_parent_directory(config.data_output_file.c_str()));
    }

    uint64_t shots_remaining = local_shots;
    uint64_t batchno = world_rank;
    while (shots_remaining) {
        const uint64_t shots_this_batch = shots_remaining < G_SHOTS_PER_BATCH ? shots_remaining : G_SHOTS_PER_BATCH;
        run_batch(program, shots_this_batch);
        write_stats(batchno);

        is_recording_stim_instructions = false;

        shots_remaining -= shots_this_batch;
        batchno += world_size;
    }
    shot_histogram = histogram_reduce(shot_histogram);
    return shot_histogram;
}

inline void
FullSystemSimulator::load_subroutine(std::string name, const qes::Program<>& program) {
    subroutine_map[name] = program;
}

inline void
FullSystemSimulator::execute_routine(const qes::Program<>& program) {
    program_status_t status(current_shots);
    while (status.pc < program.size()) {
        const qes::Instruction<>& instruction = program.at(status.pc);
        // Check if the instruction is requesting the execution of the
        // microcode.
        if (instruction.get_name() == "call") {
            // Switch into the subroutine.
            std::string subroutine_name = instruction.get<std::string>(0);
            execute_routine(subroutine_map.at(subroutine_name));
            status.pc++;
        } else {
            read_next_instruction(program, status);
        }
    }
}

inline stim::simd_bits_range_ref<SIMD_WIDTH>
FullSystemSimulator::get_register(std::string r) {
    // Handle certain registers differently. If it is a special register, we may need to
    // update the value here.
    if (r == "$reoz") {
        // Set $reoz to the inverse of the bitwise OR of the 
        // last config.sse.rreoz_track_n_det events.
        size_t k = get_register_index("$reoz");
        register_file[k].clear();
        for (uint64_t e : sse.rreoz_events) {
            register_file[k] |= syndrome_table[e];
        }
        register_file[k].invert_bits();
    }
    return register_file[get_register_index(r)];
}

inline void
FullSystemSimulator::recalibrate_timing() {
    // Compute min time when taking account delta.
    fp_t min_time_delta = 0.0;
    for (const auto& p : shot_time_delta_map) {
        min_time_delta = std::min(min_time_delta, p.second);
    }
    if (min_time_delta >= 0.0) return;
    // Otherwise, change everything up.
    std::map<uint64_t, fp_t> new_time_delta_map;
    for (uint64_t t = 0; t < current_shots; t++) {
        if (shot_time_delta_map.count(t)) {
            if (min_time_delta != shot_time_delta_map.at(t)) {
                new_time_delta_map[t] = shot_time_delta_map.at(t) - min_time_delta;
            }
        } else {
            new_time_delta_map[t] = -min_time_delta;
        }
    }
    elapsed_time += min_time_delta;
    shot_time_delta_map = std::move(new_time_delta_map);
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
