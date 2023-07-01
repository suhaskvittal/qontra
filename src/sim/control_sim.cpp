/*
 *  author: Suhas Vittal
 *  date:   20 June 2023
 * */

#include "sim/control_sim.h"

namespace qontra {

using namespace experiments;

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
    pauli_frames(128, G_SHOTS_PER_BATCH),
    event_history(4096, G_SHOTS_PER_BATCH),
    obs_buffer(128, G_SHOTS_PER_BATCH),
    obs_buffer_max_written(0),
    // Microarchitecture
    program(),
    pc(G_SHOTS_PER_BATCH, 64),
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
    rng(0),
    flag_canonical_circuit(false),
    canonical_circuit(),
    is_fast_forwarding(false),
    apply_pending_errors(false)
{
    program.fill({"nop", {}, {}});
    for (uint i = 0; i < prog.size(); i++) {
        program[i] = prog[i];
    }
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

    rng.seed(G_BASE_SEED + world_rank);

    trial_done.clear();
    uint64_t shots_left = shots;
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
            std::cout << "[ BATCH shots = "
                << shots_in_curr_batch << " ]\n";
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
                __n_trials_killed += trial_done.popcnt();
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
        shots_left -= shots_in_curr_batch;
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

    latency.clear();

    pc.clear();
    
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
ControlSimulator::build_canonical_circuit() {
    // Build canonical circuit by running through the simulator
    // once through the program.
    flag_canonical_circuit = true;
    this->run(1);
    flag_canonical_circuit = false;
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
        if (params.enable_fast_forwarding && !flag_canonical_circuit) {
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
        if (inst.exclude_trials.count(t)) {
            id_qex_valid[t] = 0;
        }
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
    stim::simd_bit_table events_t = event_history.transposed();
    stim::simd_bit_table obs_t = obs_buffer.transposed();

    // As we can execute quantum instructions in a batch, track
    // which instructions are being executed in which trials.
    std::map<Instruction, std::set<uint>>   instruction_to_trials;

    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if (qex_stall[t] || !id_qex_valid[t])   continue;
        auto inst = program[id_pc[t].u64[0]];

        if (IS_NOP_LIKE.count(inst.name))    continue;

        // If the instruction requires interacting with any qubits,
        // check if any operands are busy. If so, stall the pipeline.
        bool stall_pipeline = false;
        if (HAS_QUBIT_OPERANDS.count(inst.name)) {
            for (uint x : inst.operands) {
                if (qex_qubit_busy[t].u64[x] > 0) {
                    qex_stall[t] = 1;
                    qex_rt_valid[t] = 0;
                    stall_pipeline = true;
                    break;
                }
            }
        } else if (IS_FENCE.count(inst.name)) {
            for (uint x = 0; x < n_qubits; x++) {
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
        uint64_t    br_pc;
        bool        br_taken = false;

        if (inst.name == "decode" && decoder != nullptr) {
            syndrome_t s(events_t[t]);
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
        } else if (inst.name == "jmp") {
            br_pc = inst.operands[0];
            br_taken = true;
        } else if (inst.name == "brdb") {
            br_pc = inst.operands[0];
            br_taken = decoder_busy[t].not_zero();
        } else if (inst.name == "braspc") {
            br_pc = inst.operands[0];
        } else if (inst.name == "brospc") {
            br_pc = inst.operands[0];
        } else if (inst.name == "dfence") {
            if (decoder_busy[t].not_zero()) {
                qex_stall[t] = 1;
                qex_rt_valid[t] = 0;
            }
        } else if (inst.name == "savem") {
            vlw_t x = to_vlw(obs_t[t], obs_buffer_max_written);
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

        if (br_taken) {
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

        if (inst.name == "h") {
            qsim->H(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("H", inst.operands);
            }
        } else if (inst.name == "x") {
            qsim->X(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("X", inst.operands);
            }
        } else if (inst.name == "z") {
            qsim->Z(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("Z", inst.operands);
            }
        } else if (inst.name == "s") {
            qsim->S(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("S", inst.operands);
            }
        } else if (inst.name == "cx") {
            qsim->CX(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("CX", inst.operands);
            }
        } else if (inst.name == "mrc" || inst.name == "mnrc") {
            qsim->M(inst.operands, inst.name == "mrc");
            if (flag_canonical_circuit && inst.name == "mrc") {
                canonical_circuit.append_op("M", inst.operands);
            }
        } else if (inst.name == "reset") {
            qsim->R(inst.operands);
            if (flag_canonical_circuit) {
                canonical_circuit.append_op("R", inst.operands);
            }
        } else if (inst.name == "event") {
            const uint k = inst.operands[0];
            event_history[k].clear();
            const uint len = qsim->get_record_size();
            for (uint i = 1; i < inst.operands.size(); i++) {
                event_history[k] ^= qsim->record_table[len - inst.operands[i]];
            }

            if (flag_canonical_circuit) {
                std::vector<uint> stim_operands;
                for (uint i = 1; i < inst.operands.size(); i++) {
                    stim_operands.push_back(inst.operands[i]
                                                | stim::TARGET_RECORD_BIT);
                }
                canonical_circuit.append_op("DETECTOR", stim_operands);
            }
        } else if (inst.name == "obs") {
            const uint k = inst.operands[0];
            obs_buffer[k].clear();
            const uint len = qsim->get_record_size();
            for (uint i = 1; i < inst.operands.size(); i++) {
                obs_buffer[k] ^= qsim->record_table[len - inst.operands[i]];
            }
            if (k+1 > obs_buffer_max_written) obs_buffer_max_written = k+1;

            if (flag_canonical_circuit) {
                std::vector<uint> stim_operands;
                for (uint i = 1; i < inst.operands.size(); i++) {
                    stim_operands.push_back(inst.operands[i]
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
                    event_history[k][t] = 0;
                } else if (inst.name == "obs") {
                    const uint k = inst.operands[0];
                    obs_buffer[k][t] = 0;
                } else if (inst.name == "xorfr") {
                    const uint k1 = inst.operands[0];
                    const uint k2 = inst.operands[1];
                    obs_buffer[k2][t] ^= pauli_frames[k1][t];
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
                if (key == "mrc" || key == "mnrc")  key = "m";
                if (!params.timing.op1q.count(inst.name))   continue;
                for (uint x : inst.operands) {
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

            if (flag_canonical_circuit) {
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
        if (inst.name == "mrc" || inst.name == "mnrc")  key = "m";
        if (!params.errors.op1q.count(key)) return;
        std::vector<fp_t> e;
        for (uint x : inst.operands) {
            e.push_back(params.errors.op1q[key][x]);
            if (flag_canonical_circuit) {
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
            qsim->eX(inst.operands, e);
        } else {
            qsim->eDP1(inst.operands, e);
        }
    }
}

void
ControlSimulator::apply_periodic_error(fp_t t) {
    std::vector<uint> q;
    if (params.simulate_periodic_as_dpo_and_dph 
        && !flag_canonical_circuit) 
    {
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

            if (flag_canonical_circuit) {
                canonical_circuit.append_op("DEPOLARIZE1", {i}, x);
            }
        }
        qsim->eDP1(q, e);
    }
}

}   // qontra
