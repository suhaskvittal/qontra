/*
 *  author: Suhas Vittal
 *  date:   20 June 2023
 * */

#include "sim/control_sim.h"

namespace qontra {

using namespace experiments;

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
        decoder_busy.clear();

        csim.reset_sim();
        csim.shots = shots_in_curr_batch;
        pauli_frames.clear();
        obs_buffer.clear();
        obs_buffer_size.clear();

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
                std::cout << "No. of trials done: " << 
                    (shots_in_curr_batch - trial_done.popcnt()) << "\n";
            }
            RT();
            QEX();
            ID();
            IF();

            // Update timing.
            uint64_t tt = (uint64_t)((1.0/params.clock_frequency) * 1e9);
            for (uint64_t i = 0; i < shots_in_curr_batch; i++) {
                if (!trial_done[i]) continue;
                latency[i].u64[0] += tt;

                if (decoder_busy[i].u64[0] <= tt) {
                    decoder_busy[i].u64[0] = 0;
                } else {
                    decoder_busy[i].u64[0] -= tt;;
                }
                for (uint64_t q = 0; q < n_qubits; q++) {
                    if (qex_qubit_busy[i].u64[q] <= tt) {
                        qex_qubit_busy[i].u64[q] = 0;
                    } else {
                        qex_qubit_busy[i].u64[q] -= tt;
                    }
                }
            }
            t_since_last_periodic_error -= tt;
            if (t_since_last_periodic_error <= 0) {
                apply_periodic_error();
            }
            // Check if we need to abort any trials.
            auto rt = timer.clk_end();
            if (rt > params.kill_batch_after_time_elapsed) {
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
    } else {
        n_trials_killed = __n_trials_killed;
        latency_max = __latency_max;
        latency_sum = __latency_sum;
        latency_sqr_sum = __latency_sqr_sum;
    }
    latency_mean = MEAN(latency_sum, shots);
    latency_std = STD(latency_mean, latency_sqr_sum, shots);
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

    // Keep the syndromes in the following data structure -- hopefully
    // they are cached.
    //
    // Assume these are already detection events.
    stim::simd_bit_table syndromes = csim.record_table.transposed();
    // As we can execute quantum instructions in a batch, track
    // which instructions are being executed in which trials.
    std::map<Instruction, std::set<uint>>   instruction_to_trials;

    for (uint64_t t = 0; t < shots_in_curr_batch; t++) {
        if (qex_stall[t] || !id_qex_valid[t])   continue;
        auto inst = program[id_pc[t].u64[0]];
        // Check if any operands are busy. If so, stall the pipeline.
        bool stall_pipeline = false;
        for (uint x : inst.operands) {
            if (qex_qubit_busy[t].u64[x] > 0) {
                qex_stall[t] = 1;
                qex_rt_valid[t] = 0;
                stall_pipeline = true;
            }
        }
        if (stall_pipeline) continue;
        // These instructions operate on a trial by trial basis.
        if (inst.name == "decode") {
            syndrome_t s(syndromes[t]);
            auto res = decoder->decode_error(s);
            uint64_t exec_time = (uint64_t)res.exec_time;
            // Decoder may get backlogged if it does not complete on time.
            //
            // But this is fine.
            decoder_busy[t].u64[0] += exec_time;
            // Update Pauli frames
            for (uint i = 0; i < res.corr.size(); i++) {
                pauli_frames[t][i] ^= res.corr[i];
            }
        } else if (inst.name == "brdb") {
            if (decoder_busy[t].not_zero()) {
                // Set PC and invalidate ID.
                id_qex_valid[t] = 0;
                pc[t].u64[0] = inst.operands[0];
            }
        } else if (inst.name == "dfence") {
            if (decoder_busy[t].not_zero()) {
                id_stall[t] = 1;
                qex_rt_valid[t] = 0;
            }
        } else if (inst.name == "event") {
            uint src = inst.operands[0];
            uint dst = inst.operands[1];
            csim.xor_record_with(src, dst);
        } else if (inst.name == "obs") {
            uint8_t obs = 0;
            uint len = csim.get_record_size();
            for (uint offset : inst.operands) {
                obs ^= csim.record_table[len-offset][t];
            }
            obs_buffer[t][obs_buffer_size[t].u64[0]++] = obs;
        } else if (inst.name == "savem") {
            vlw_t x = to_vlw(obs_buffer[t], obs_buffer_size[t].u64[0]);
            prob_histograms[x]++;
            obs_buffer[t].clear();
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
                if (!pair.second.count(t))  std::cout << " " << t;
            }
            std::cout << "\n";
        }
        // Create snapshot of Clifford simulator state.
        csim.snapshot();
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
            csim.H(inst.operands);
        } else if (inst.name == "x") {
            csim.X(inst.operands);
        } else if (inst.name == "z") {
            csim.Z(inst.operands);
        } else if (inst.name == "s") {
            csim.S(inst.operands);
        } else if (inst.name == "cx") {
            csim.CX(inst.operands);
        } else if (inst.name == "mrc" || inst.name == "mnrc") {
            csim.M(inst.operands, inst.name == "mrc");
        } else if (inst.name == "reset") {
            csim.R(inst.operands);
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
                csim.rollback_at_trial(t);
            } else {
                if (inst.name == "CX") {
                    for (uint i = 0; i < inst.operands.size(); i += 2) {
                        uint x = inst.operands[i];
                        uint y = inst.operands[i+1];
                        fp_t tt = params.timing.op2q[inst.name][std::make_pair(x, y)];
                        qex_qubit_busy[t].u64[x] = tt;
                        qex_qubit_busy[t].u64[y] = tt;
                    }
                } else {
                    if (params.timing.op1q.count(inst.name)) {
                        for (uint x : inst.operands) {
                            fp_t tt = params.timing.op1q[inst.name][x];
                            qex_qubit_busy[t].u64[x] = tt;
                        }
                    } else if (inst.name == "mrc" || inst.name == "mnrc") {
                        for (uint x : inst.operands) {
                            fp_t tt = params.timing.op1q["m"][x];
                            qex_qubit_busy[t].u64[x] = tt;
                        }
                    }
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
    if (inst.name == "CX") {
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
        }
        csim.eLT(inst.operands, lt);
        csim.eLI(inst.operands, li);
        csim.eDP2(inst.operands, dp2);
    } else {
        std::string key = inst.name;
        if (inst.name == "mrc" || inst.name == "mnrc")  key = "m";
        if (!params.errors.op1q.count(key)) return;
        std::vector<fp_t> e;
        for (uint x : inst.operands) {
            e.push_back(params.errors.op1q[key][x]);
        }
        if (key == "m" || key == "reset") {
            csim.eX(inst.operands, e);
        } else {
            csim.eDP1(inst.operands, e);
        }
    }
}

void
ControlSimulator::apply_periodic_error() {
    const fp_t t = params.apply_periodic_errors_at_t;
    std::vector<uint> q;
    if (params.simulate_periodic_as_dpo_and_dph) {
        std::vector<fp_t> e1, e2;
        for (uint i = 0; i < n_qubits; i++) {
            q.push_back(i);
            fp_t x = 1 - exp(-t / params.timing.t1[i]);
            fp_t y = 1 - exp(-t / params.timing.t2[i]);
            e1.push_back(x);
            e2.push_back(y);
        }
        csim.eDPO(q, e1);
        csim.eDPH(q, e2);
    } else {
        std::vector<fp_t> e;
        for (uint i = 0; i < n_qubits; i++) {
            fp_t mt = 0.5 * (params.timing.t1[i] + params.timing.t2[i]);
            q.push_back(i);
            e.push_back(1 - exp(-t / mt));
        }
        csim.eDP1(q, e);
    }
}

}   // qontra
