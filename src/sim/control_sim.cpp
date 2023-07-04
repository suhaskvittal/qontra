/*
 *  author: Suhas Vittal
 *  date:   20 June 2023
 * */

#include "sim/control_sim.h"

namespace qontra {

namespace ctrlsim {

uint64_t    EVENT_HISTORY_SIZE = 1024*8;
uint64_t    OBS_BUFFER_SIZE = 128;
uint64_t    TRACES_PER_SHOT = 1;

}

using namespace experiments;
using namespace ctrlsim;

#define N_TRACES (G_SHOTS_PER_BATCH*TRACES_PER_SHOT)

ControlSimulator::ControlSimulator(
        uint n_qubits, 
        const schedule_t& prog,
        StateSimulator* qsim)
    // Stats
    :latency(G_SHOTS_PER_BATCH, 64),
    sim_time(0),
    // Simulation state tracking
    decoder_busy(G_SHOTS_PER_BATCH, 64),
    trial_done(G_SHOTS_PER_BATCH),
    // Simulation structures
    decoder(nullptr),
    qsim(qsim),
    pauli_frames(OBS_BUFFER_SIZE, G_SHOTS_PER_BATCH),
    event_history(EVENT_HISTORY_SIZE, G_SHOTS_PER_BATCH),
    obs_buffer(OBS_BUFFER_SIZE, G_SHOTS_PER_BATCH),
    obs_buffer_max_written(0),
    event_trace(N_TRACES, EVENT_HISTORY_SIZE+OBS_BUFFER_SIZE),
    event_trace_index(G_SHOTS_PER_BATCH, 64),
    event_trace_max_written(0),
    // Microarchitecture
    program(),
    pc(G_SHOTS_PER_BATCH, 64),
    qubit_locks(G_SHOTS_PER_BATCH, n_qubits),
    // Selective Syndrome Extraction
    sig_m_spec(4096, G_SHOTS_PER_BATCH),
    val_m_spec(4096, G_SHOTS_PER_BATCH),
    // IF io
    if_stall(G_SHOTS_PER_BATCH),
    if_pc(G_SHOTS_PER_BATCH, 64),
    if_id_valid(G_SHOTS_PER_BATCH),
    // ID io
    id_stall(G_SHOTS_PER_BATCH),
    id_pc(G_SHOTS_PER_BATCH, 64),
    id_qex_valid(G_SHOTS_PER_BATCH),
    // QEX io
    qex_stall(G_SHOTS_PER_BATCH),
    qex_rt_valid(G_SHOTS_PER_BATCH),
    qex_pc(G_SHOTS_PER_BATCH, 64),
    qex_qubit_busy(G_SHOTS_PER_BATCH, 64*n_qubits),
    // RT io
    rt_stall(G_SHOTS_PER_BATCH),
    // Other
    n_qubits(n_qubits),
    is_fast_forwarding(false),
    apply_pending_errors(false)
{
    program.fill({"nop", {}, {}});
    for (uint i = 0; i < prog.size(); i++) {
        program[i] = prog[i];
    }
}

void
ControlSimulator::build_error_model() {
    measurement_order_to_loc.clear();
    measurement_count = 0;

    bool use_mpi = G_USE_MPI;
    G_USE_MPI = false;

    is_building_canonical_circuit = true;
    run(1);
    is_building_canonical_circuit = false;

    G_USE_MPI = use_mpi;
}


void
ControlSimulator::run(uint64_t shots) {
    int world_rank = 0, world_size = 1;
    uint64_t local_shots = shots;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        local_shots = shots / world_size;
        if (world_rank == 0) local_shots += shots % world_size; 
    }

    uint64_t __n_trials_killed = 0;
    uint64_t __latency_max = 0;
    uint64_t __sim_time = 0;
    uint64_t latency_sum, __latency_sum = 0;
    uint64_t latency_sqr_sum, __latency_sqr_sum = 0;

    qsim->set_seed(G_BASE_SEED + world_rank);

    trial_done.clear();

    uint64_t shots_left = local_shots;
    uint bno = world_rank;   // Batch number

    if (params.save_syndromes_to_file && world_rank == 0) {
        // Delete the folder and create it again to create a fresh trace
        // folder.
        std::filesystem::path folder_path(params.syndrome_output_folder);
        std::filesystem::remove_all(folder_path);
        safe_create_directory(folder_path);
    }

    while (shots_left) {
        shots_in_curr_batch = shots_left < G_SHOTS_PER_BATCH 
                                ? shots_left : G_SHOTS_PER_BATCH;
        if (shots_in_curr_batch < G_SHOTS_PER_BATCH) {
            for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
                trial_done[t] = 1;
            }
        } else {
            trial_done.invert_bits();
        }
        // Reset structures.
        clear();
        // Start simulation
        timer.clk_start();

        if (params.verbose) {
            std::cout << "[ BATCH shots = " << shots_in_curr_batch << " ]\n";
        }
        
        int64_t t_since_last_periodic_error = params.apply_periodic_errors_at_t;
        while (trial_done.not_zero()) {
            if (params.verbose) {
                std::cout << "\tPC(t = 0, T = "
                    << latency[0].u64[0] << ") = " 
                    << pc[0].u64[0] << "\t"
                    << "@ " << program[pc[0].u64[0]].str() << "\n";
                std::cout << "\tIs fast-forwarding: " << is_fast_forwarding << "\n";
            }
            RT();
            QEX();
            ID();
            IF();

            // Update timing.
            uint64_t tt = (uint64_t)((1.0/params.clock_frequency) * 1e9);
            for (uint64_t i = 0; i < shots_in_curr_batch; i++) {
                if (!trial_done[i]) continue;
                // Do not count timing statistics due to fastforwarding.
                if (!is_fast_forwarding) {
                    latency[i].u64[0] += tt;
                }

                if (decoder_busy[i].u64[0] <= tt || is_fast_forwarding) {
                    decoder_busy[i].u64[0] = 0;
                } else {
                    decoder_busy[i].u64[0] -= tt;;
                }
                for (uint64_t q = 0; q < n_qubits; q++) {
                    if (qex_qubit_busy[i].u64[q] <= tt
                        || is_fast_forwarding) 
                    {
                        qex_qubit_busy[i].u64[q] = 0;
                    } else {
                        qex_qubit_busy[i].u64[q] -= tt;
                    }
                }
            }

            // Apply any periodic errors if necessary.
            t_since_last_periodic_error -= tt;
            if (t_since_last_periodic_error <= 0) {
                if (is_fast_forwarding) {
                    apply_periodic_error(
                            params.ff_periodic_error_assume_time_elapsed);
                    t_since_last_periodic_error = 
                            params.ff_apply_periodic_errors_at_t;
                } else {
                    apply_periodic_error(params.apply_periodic_errors_at_t);
                    t_since_last_periodic_error = 
                            params.apply_periodic_errors_at_t;
                }
            } else if (apply_pending_errors) {
                // This is a transitionary period between FF and non-FF (or
                // vice versa).
                if (t_since_last_periodic_error < 0) {
                    t_since_last_periodic_error = 0;
                }
                uint64_t delta;
                if (is_fast_forwarding) {
                    delta = params.apply_periodic_errors_at_t
                                - t_since_last_periodic_error;
                    delta *= params.ff_periodic_error_assume_time_elapsed
                                / params.ff_apply_periodic_errors_at_t;
                    t_since_last_periodic_error = 
                            params.apply_periodic_errors_at_t;
                } else {
                    delta = params.ff_periodic_error_assume_time_elapsed 
                                - t_since_last_periodic_error;
                    t_since_last_periodic_error = 
                            params.ff_apply_periodic_errors_at_t;
                }
                apply_periodic_error(delta);
            }
            apply_pending_errors = false;

            // Check if we need to abort any trials.
            auto rt = timer.clk_end();
            if (params.verbose) {
                std::cout << "\tTime elapsed since batch start: " 
                    << rt*1e-9 << "\n";
            }
            if (rt*1e-9 > params.kill_batch_after_time_elapsed) {
                std::cout << "Killed " << trial_done.popcnt()
                        << " trials because walltime exceeded "
                        << rt*1e-9 << "s.\n";
                __n_trials_killed += trial_done.popcnt();
                std::cout << "[ ALERT ] Killed "
                        << trial_done.popcnt() << 
                        " since walltime of " 
                        << rt*1e-9
                        << "s has been exceeded.\n";
                trial_done.clear();
            }
        }
        // Update stats.
        for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
            uint64_t x = latency[t].u64[0];
            __latency_sum += x;
            __latency_sqr_sum += SQR(x);
            if (x > __latency_max)  __latency_max = x;
        }
        auto rt = timer.clk_end();
        __sim_time += rt;
        // Save syndrome trace if requested.
        if (params.save_syndromes_to_file) {
            std::string filename = "batch_" + std::to_string(bno) + ".dets";
            std::string output_path = params.syndrome_output_folder 
                                    + "/" + filename;

            FILE* fout = fopen(output_path.c_str(), "w");
            const uint trace_width = event_trace_max_written
                                    + params.save_observables.size();
            stim::simd_bits ref(trace_width);   // Because Stim
                                                // decides we need
                                                // this for some
                                                // reason.
            stim::write_table_data(fout, 
                            shots_in_curr_batch*TRACES_PER_SHOT,
                            trace_width,
                            ref,
                            event_trace,
                            stim::SampleFormat::SAMPLE_FORMAT_DETS,
                            'D',
                            'L',
                            event_trace_max_written);
            fclose(fout);
        }
        shots_left -= shots_in_curr_batch;
        bno += world_size;
    }

    if (G_USE_MPI) {
        MPI_Allreduce(&__n_trials_killed, &n_trials_killed, 1, MPI_UNSIGNED_LONG,
                    MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&__latency_max, &latency_max, 1, MPI_UNSIGNED_LONG,
                    MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&__latency_sum, &latency_sum, 1, MPI_UNSIGNED_LONG,
                    MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&__latency_sqr_sum, &latency_sqr_sum, 1, MPI_UNSIGNED_LONG,
                    MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&__sim_time, &sim_time, 1, MPI_UNSIGNED_LONG,
                    MPI_MAX, MPI_COMM_WORLD);
        // Accumulate probability histograms. This will be a bit hard.
        std::map<vlw_t, uint64_t> local_ph(prob_histograms);
        for (uint r = 0; r < world_size; r++) {
            if (world_rank == r) {
                for (auto pair : local_ph) {
                    // Transmit which vlw_t is being sent.
                    vlw_t x = pair.first;
                    const uint xw = x.size();
                    const uint64_t count = pair.second;
                    for (uint s = 0; s < world_size; s++) {
                        if (s == r) continue;
                        // Transmit size of vlw_t first.
                        MPI_Bsend(&xw, 1, MPI_UNSIGNED, s, 0, MPI_COMM_WORLD);
                        // Now, transmit the vlw_t.
                        MPI_Bsend(&x[0], xw, MPI_UNSIGNED_LONG, s, 1, MPI_COMM_WORLD);
                        // Finally, transmit the count.
                        MPI_Bsend(&count, 1, MPI_UNSIGNED_LONG, s, 2, MPI_COMM_WORLD);
                    }
                }
                for (uint s = 0; s < world_size; s++) {
                    if (s == r) continue;
                    // Finally, transmit that we are done.
                    uint done = 0;
                    MPI_Bsend(&done, 1, MPI_UNSIGNED, s, 0, MPI_COMM_WORLD);
                }
            } else {
                while (true) {
                    uint xw;
                    MPI_Recv(&xw, 1, MPI_UNSIGNED, r, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (!xw)    break;
                    vlw_t x(xw);
                    MPI_Recv(&x[0], xw, MPI_UNSIGNED_LONG, r, 1,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    uint64_t count;
                    MPI_Recv(&count, 1, MPI_UNSIGNED_LONG, r, 2,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (!prob_histograms.count(x))  prob_histograms[x] = 0;
                    prob_histograms[x] += count;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    } else {
        n_trials_killed = __n_trials_killed;
        latency_max = __latency_max;
        latency_sum = __latency_sum;
        latency_sqr_sum = __latency_sqr_sum;
        sim_time = __sim_time;
    }
    latency_mean = MEAN(latency_sum, shots);
    latency_std = STD(latency_mean, latency_sqr_sum, shots);
}

void
ControlSimulator::clear() {
    decoder_busy.clear();

    qsim->reset_sim();
    qsim->shots = shots_in_curr_batch;
    pauli_frames.clear();
    event_history.clear();
    obs_buffer.clear();
    obs_buffer_max_written = 0;

    event_trace.clear();
    event_trace_index.clear();
    event_trace_max_written = 0;

    latency.clear();

    pc.clear();
    qubit_locks.clear();
    
    if_stall.clear();
    if_id_valid.clear();
    if_pc.clear();
    id_stall.clear();
    id_qex_valid.clear();
    id_pc.clear();
    qex_stall.clear();
    qex_rt_valid.clear();
    qex_pc.clear();
    qex_qubit_busy.clear();
    rt_stall.clear();
}

void
ControlSimulator::IF() {
    if (params.verbose) {
        std::cout << "\t[ IF ]\n";
    }

    if_stall = id_stall;
    // Set all to 1.
    if_id_valid.clear();
    if_id_valid.invert_bits();
    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if (!trial_done[t]) continue;
        if (if_stall[t])    continue;
        // Otherwise, update the PC and do any IO.
        if_pc[t].u64[0] = pc[t].u64[0];

        // Check if this is a fastforwarding instruction.
        auto& inst = program[pc[t].u64[0]];
        if (params.enable_fast_forwarding && !is_building_canonical_circuit) {
            if (inst.name == "ffstart") {
                // Apply any pending periodic errors.
                apply_pending_errors = true;
                is_fast_forwarding = true;
            } else if (inst.name == "ffend") {
                apply_pending_errors = true;
                is_fast_forwarding = false;
            }
        }
        pc[t].u64[0]++;
    }
}

void
ControlSimulator::ID() {
    if (params.verbose) {
        std::cout << "\t[ ID ]\n";
    }

    id_stall = qex_stall;
    id_qex_valid = if_id_valid;
    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if (!trial_done[t]) continue;
        if (id_stall[t] || !if_id_valid[t])   continue;
        id_pc[t].u64[0] = if_pc[t].u64[0];
        // Check if the instruction is excluded from this trial.
        // If so, invalidate the pipeline latch.
        auto inst = program[id_pc[t].u64[0]];
    }
}

void
ControlSimulator::QEX() {
    if (params.verbose) {
        std::cout << "\t[ QEX ]\n";
    }

    qex_stall = rt_stall;
    qex_rt_valid = id_qex_valid;

    // Keep the following data structures' transposes here.
    // The assumption is that by temporal locality, we will be
    // accessing the contents of these structures each trial (assuming
    // that most trials are at the same point in the program). So,
    // this should maximize cache accesses.
    stim::simd_bit_table events_tp = event_history.transposed();
    stim::simd_bit_table obs_tp = obs_buffer.transposed();

    // As we can execute quantum instructions in a batch, track
    // which instructions are being executed in which trials.
    std::map<Instruction, std::set<uint>>   instruction_to_trials;

    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if (qex_stall[t] || !id_qex_valid[t])   continue;
        auto inst = program[id_pc[t].u64[0]];

        if (IS_NOP_LIKE.count(inst.name))    continue;

        auto qubit_operands = inst.get_qubit_operands();
        // If the instruction requires interacting with any qubits,
        // check if any operands are busy. If so, stall the pipeline.
        bool stall_pipeline = false;
        if (is_building_canonical_circuit && IS_LOCKING.count(inst.name)) {
            continue;
        } else if (IS_FENCE.count(inst.name)) {
            for (uint x = 0; x < n_qubits; x++) {
                if (qex_qubit_busy[t].u64[x] > 0) {
                    qex_stall[t] = 1;
                    qex_rt_valid[t] = 0;
                    stall_pipeline = true;
                    break;
                }
            }
        } else {
            for (auto x : qubit_operands) {
                if (qex_qubit_busy[t].u64[x] > 0) {
                    qex_stall[t] = 1;
                    qex_rt_valid[t] = 0;
                    stall_pipeline = true;
                    break;
                }
            }
        }
        if (stall_pipeline) continue;

        // These instructions operate on a trial by trial basis.
        // 
        // See instruction.h for details about what these instructions
        // do.
        uint64_t    br_pc;
        bool        br_taken = false;

        if (inst.name == "decode") {
            syndrome_t s(events_tp[t]);
            if (decoder != nullptr) {
                auto res = decoder->decode_error(s);
                uint64_t exec_time = (uint64_t)res.exec_time;
                // Decoder may get backlogged if it does not complete on time.
                //
                // But this is fine.
                if (!params.decoder_is_ideal) {
                    decoder_busy[t].u64[0] += exec_time;
                }
                // Update Pauli frames
                uint offset = inst.operands[0];
                for (uint i = 0; i < res.corr.size(); i++) {
                    pauli_frames[i+offset][t] ^= res.corr[i];
                }
            }

            if (params.save_syndromes_to_file) {
                // Write syndrome to event_trace.
                uint64_t& trace_index = event_trace_index[t].u64[0];
                uint64_t i = t*TRACES_PER_SHOT + trace_index;
                event_trace[i].clear();
                event_trace[i] |= s;
                for (uint j = 0; j < params.save_observables.size(); j++) {
                    uint obs = params.save_observables[j];
                    event_trace[i][event_trace_max_written+j] = obs_tp[t][j];
                }
                trace_index++;
            }
        } else if (inst.name == "lockq") {
            for (uint i : inst.operands) {
                qubit_locks[t][i] = 1;
            }
        } else if (inst.name == "unlockq") {
            for (uint i : inst.operands) {
                qubit_locks[t][i] = 0;
            }
        } else if (inst.name == "jmp") {
            br_pc = inst.operands[0];
            br_taken = true;
        } else if (inst.name == "brdb") {
            br_pc = inst.operands[0];
            br_taken = decoder_busy[t].not_zero();
        } else if (inst.name == "brifmspc") {
            br_pc = inst.operands[0];
            bool all_speculated = true;
            for (uint i = 1; i < inst.operands.size(); i++) {
                uint j = inst.operands[i];
                all_speculated &= sig_m_spec[j][t];
                // Premptively place the results.
                qsim->record_table[j][t] = val_m_spec[j][t];
            }
            br_taken = all_speculated;
        } else if (inst.name == "lockqifmspc") {
            bool all_speculated = true;
            for (uint i = 1; i < inst.operands.size(); i++) {
                uint j = inst.operands[i];
                all_speculated &= sig_m_spec[j][t];
            }
            // Lock the qubit operand.
            qubit_locks[t][inst.operands[0]] = all_speculated;
        } else if (inst.name == "dfence") {
            if (decoder_busy[t].not_zero()) {
                qex_stall[t] = 1;
                qex_rt_valid[t] = 0;
            }
        } else if (inst.name == "savem") {
            vlw_t x = to_vlw(obs_tp[t], obs_buffer_max_written);
            prob_histograms[x]++;
        } else if (inst.name == "done") {
            trial_done[t] = 0;
        } else {
            // These should be instructions that benefit from operating
            // on multiple trials at once (i.e. quantum gates).
            //
            // We can execute the gates for all trials, and then rollback
            // the changes for any trials that did not want to execute the
            // gate.
            instruction_to_trials[inst].insert(t);
        }

        if (br_taken && !is_building_canonical_circuit) {
            if_id_valid[t] = 0;
            id_qex_valid[t] = 0;
            pc[t].u64[0] = br_pc;
        }
    }
    // Now execute the instructions in instruction_to_trials.
    if (params.verbose) {
        std::cout << "\t\tInstructions executed in batch:\n";
    }
    for (auto pair : instruction_to_trials) {
        auto inst = pair.first;
        if (params.verbose) {
            std::cout << "\t\t\t" << inst.str() << "\tT !=";
            for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
                if (!pair.second.count(t)) {
                    std::cout << " " << t;
                }
            }
            std::cout << "\n";
        }
        // Create snapshot of Clifford simulator state.
        qsim->snapshot();
        // Execute the operation, alongside any errors.
        //
        // Before operation errors:
        const std::set<std::string> apply_error_before{
            "mrc",
            "mnrc"
        };

        if (apply_error_before.count(inst.name)) {
            apply_gate_error(inst);
        }

        stim::simd_bit_table    event_history_cpy(event_history);
        stim::simd_bit_table    sig_m_spec_cpy(sig_m_spec);
        stim::simd_bit_table    val_m_spec_cpy(val_m_spec);

        if (inst.name == "h") {
            qsim->H(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("H", inst.operands);
            }
        } else if (inst.name == "x") {
            qsim->X(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("X", inst.operands);
            }
        } else if (inst.name == "z") {
            qsim->Z(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("Z", inst.operands);
            }
        } else if (inst.name == "s") {
            qsim->S(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("S", inst.operands);
            }
        } else if (inst.name == "cx") {
            qsim->CX(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("CX", inst.operands);
            }
        } else if (inst.name == "mrc") {
            auto qubits = inst.get_qubit_operands();
            qsim->M(qubits, inst.operands[0]);
            if (is_building_canonical_circuit) {
                for (uint i = 0; i < inst.operands.size()-1; i++) {
                    measurement_order_to_loc[inst.operands[0]+i] = 
                        measurement_count++;
                }
                canonical_circuit.append_op("M", qubits);
            }
        } else if (inst.name == "mnrc") {
            qsim->M(inst.operands, -1);
        } else if (inst.name == "reset") {
            qsim->R(inst.operands);
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("R", inst.operands);
            }
        } else if (inst.name == "event") {
            const uint k = inst.operands[0];
            event_history[k].clear();
            for (uint i = 1; i < inst.operands.size(); i++) {
                uint j = inst.operands[i];
                event_history[k] ^= qsim->record_table[j];
            }
            if (k+1 > event_trace_max_written) event_trace_max_written = k+1;

            if (is_building_canonical_circuit) {
                std::vector<uint> stim_operands;
                for (uint i = 1; i < inst.operands.size(); i++) {
                    uint j = measurement_order_to_loc[inst.operands[i]];
                    stim_operands.push_back((measurement_count - j)
                                                | stim::TARGET_RECORD_BIT);
                }
                canonical_circuit.append_op("DETECTOR", stim_operands);
                continue;
            }
            // Do not execute the below if we are building a canonical
            // circuit.
            if (params.speculate_measurements && inst.operands.size() == 3) {
                uint m1 = inst.operands[1];
                uint m2 = inst.operands[2];
                if (m1 > m2)    std::swap(m1, m2);
                // We are assuming that a future event will have the
                // same "stride" as this event.
                uint delta = m2 - m1;
                uint m3 = m2 + delta;

                // If the event is 0, then speculate that the next
                // measurement will remain the same.
                sig_m_spec[m3].clear();
                val_m_spec[m3].clear();

                sig_m_spec[m3] |= event_history[k];
                sig_m_spec[m3].invert_bits();

                val_m_spec[m3] |= qsim->record_table[m2];
            }
        } else if (inst.name == "obs") {
            const uint k = inst.operands[0];
            obs_buffer[k].clear();
            for (uint i = 1; i < inst.operands.size(); i++) {
                obs_buffer[k] ^= qsim->record_table[inst.operands[i]];
            }
            if (k+1 > obs_buffer_max_written) obs_buffer_max_written = k+1;

            if (is_building_canonical_circuit) {
                std::vector<uint> stim_operands;
                for (uint i = 1; i < inst.operands.size(); i++) {
                    uint j = measurement_order_to_loc[inst.operands[i]];
                    stim_operands.push_back((measurement_count - j)
                                                | stim::TARGET_RECORD_BIT);
                }
                canonical_circuit.append_op("OBSERVABLE_INCLUDE", 
                        stim_operands, k);
            }
        } else if (inst.name == "xorfr") {
            const uint i = inst.operands[0];
            const uint j = inst.operands[1];
            obs_buffer[j] ^= pauli_frames[i];
        } else if (inst.name == "hshift") {
            const uint x = inst.operands[0];
            qsim->shift_record_by(x);
            // Shift event history and speculation data as well.
            for (uint i = 0; i < 4096; i++) {
                if (i < x) {
                    event_history[i].clear();
                    sig_m_spec[i].clear();
                    val_m_spec[i].clear();
                } else {
                    event_history[i].swap_with(event_history[i-x]);
                    sig_m_spec[i].swap_with(sig_m_spec[i-x]);
                    val_m_spec[i].swap_with(val_m_spec[i-x]);
                }
            }

            if (is_building_canonical_circuit) {
                std::map<uint, uint> new_mmap;
                for (auto pair : measurement_order_to_loc) {
                    if (pair.first < x) continue;
                    measurement_order_to_loc[pair.first-x] = pair.second;
                }
            }
        }

        // After operation errors:
        if (!apply_error_before.count(inst.name)) {
            apply_gate_error(inst);
        }

        // Now set the qubits as busy for X amount of time.
        //
        // Or if the trial did not want to execute the operation, rollback
        // the change.
        for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
            if (!pair.second.count(t)) {
                if (params.verbose) {
                    std::cout << "\t\tRolling back state at T = " << t << "\n";
                }
                // Undo any changes
                if (inst.name == "event") {
                    const uint k = inst.operands[0];
                    event_history[k][t] = event_history_cpy[k][t];
                    if (params.speculate_measurements 
                            && inst.operands.size() == 3) 
                    {
                        uint m1 = inst.operands[1];
                        uint m2 = inst.operands[2];
                        if (m1 > m2)    std::swap(m1, m2);
                        uint delta = m2 - m1;
                        uint m3 = m2 + delta;
                        sig_m_spec[m3][t] = sig_m_spec_cpy[m3][t];
                        val_m_spec[m3][t] = val_m_spec_cpy[m3][t];
                    }
                } else if (inst.name == "obs") {
                    const uint k = inst.operands[0];
                    obs_buffer[k][t] = 0;
                } else if (inst.name == "xorfr") {
                    const uint k1 = inst.operands[0];
                    const uint k2 = inst.operands[1];
                    obs_buffer[k2][t] ^= pauli_frames[k1][t];
                } else if (inst.name == "hshift") {
                    for (uint i = 0; i < 4096; i++) {
                        event_history[i][t] = event_history_cpy[i][t];
                    }
                    qsim->rollback_at_trial(t);
                } else {
                    qsim->rollback_at_trial(t);
                }
                continue;
            }
            // Do not do any latency updates if we are fast forwarding.
            if (is_fast_forwarding) continue;
            if (inst.name == "cx") {
                for (uint i = 0; i < inst.operands.size(); i += 2) {
                    uint x = inst.operands[i];
                    uint y = inst.operands[i+1];
                    fp_t tt = params.timing.op2q[inst.name][std::make_pair(x, y)];
                    qex_qubit_busy[t].u64[x] = tt;
                    qex_qubit_busy[t].u64[y] = tt;
                }
            } else {
                std::string key = inst.name;
                std::vector<uint> operands = inst.get_qubit_operands();
                if (key == "mrc" || key == "mnrc")  key = "m";
                if (!params.timing.op1q.count(key))   continue;
                for (uint x : operands) {
                    fp_t tt = params.timing.op1q[key][x];
                    qex_qubit_busy[t].u64[x] = tt;
                }
            }
        }
    }
}

void
ControlSimulator::RT() {
}

void
ControlSimulator::apply_gate_error(Instruction& inst) {
    if (inst.name == "cx") {
        std::vector<fp_t>   dp2;
        std::vector<fp_t>   li;
        std::vector<fp_t>   lt;
        for (uint i = 0; i < inst.operands.size(); i += 2) {
            uint x = inst.operands[i];
            uint y = inst.operands[i+1];
            auto x_y = std::make_pair(x, y);
            dp2.push_back(params.errors.op2q[inst.name][x_y]);
            li.push_back(params.errors.op2q_leakage_injection[inst.name][x_y]);
            lt.push_back(params.errors.op2q_leakage_transport[inst.name][x_y]);
            
            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("L_TRANSPORT", 
                            {x, y},
                            params.errors.op2q_leakage_transport[inst.name][x_y]);
                canonical_circuit.append_op("L_ERROR", 
                            {x, y},
                            params.errors.op2q_leakage_injection[inst.name][x_y]);
                canonical_circuit.append_op("DEPOLARIZE2", 
                            {x, y},
                            params.errors.op2q[inst.name][x_y]);
            }
        }
        qsim->eLT(inst.operands, lt);
        qsim->eLI(inst.operands, li);
        qsim->eDP2(inst.operands, dp2);
    } else {
        std::string key = inst.name;
        std::vector<uint> operands = inst.get_qubit_operands();
        if (inst.name == "mrc" || inst.name == "mnrc")  key = "m";
        if (!params.errors.op1q.count(key)) return;

        std::vector<fp_t> e;
        for (uint x : operands) {
            e.push_back(params.errors.op1q[key][x]);
            if (is_building_canonical_circuit) {
                if (key == "m" || key == "reset") {
                    canonical_circuit.append_op("X_ERROR", 
                            {x}, params.errors.op1q[key][x]);
                } else {
                    canonical_circuit.append_op("DEPOLARIZE1", 
                            {x}, params.errors.op1q[key][x]);
                }
            }
        }

        if (key == "m" || key == "reset") {
            qsim->eX(operands, e);
        } else {
            qsim->eDP1(operands, e);
        }
    }
}

void
ControlSimulator::apply_periodic_error(fp_t t) {
    std::vector<uint> q;
    if (params.simulate_periodic_as_dpo_and_dph
        && !is_building_canonical_circuit) {
        std::vector<fp_t> e1, e2;
        for (uint i = 0; i < n_qubits; i++) {
            q.push_back(i);
            fp_t x = 1 - exp(-t / params.timing.t1[i]);
            fp_t y = 1 - exp(-t / params.timing.t2[i]);
            e1.push_back(x);
            e2.push_back(y);
        }
        qsim->eDPO(q, e1);
        qsim->eDPH(q, e2);
    } else {
        std::vector<fp_t> e;
        for (uint i = 0; i < n_qubits; i++) {
            fp_t mt = 0.5 * (params.timing.t1[i] + params.timing.t2[i]);
            fp_t x = 1 - exp(-t / mt);
            q.push_back(i);
            e.push_back(x);

            if (is_building_canonical_circuit) {
                canonical_circuit.append_op("DEPOLARIZE1", {i}, x);
            }
        }
        qsim->eDP1(q, e);
    }
}

}   // qontra
