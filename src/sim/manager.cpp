/*
 *  author: Suhas Vittal
 *  date:   17 August 2023
 * */

#include "sim/manager.h"

namespace qontra {

using namespace experiments;

histogram_t
SimManager::evaluate_monte_carlo(uint64_t shots) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    const uint64_t local_shots = shots / world_size
                                    + (world_rank < (shots % world_size));
    uint64_t shots_left = local_shots;

    histogram_t local_prob_histogram;

    while (shots_left) {
        const uint64_t shots_this_batch = shots_left < G_SHOTS_PER_BATCH ?
                                            shots_left : G_SHOTS_PER_BATCH;
        // Execute batch of program.
        simulate_batch(shots_this_batch);
        // Record obs buffer. 
        uint n_words_per_obs = (max_obs_written >> 6) + 1;

        auto obt_tr = obs_buffer_table.transposed();
        for (uint64_t t = 0; t < shots_this_batch; t++) {
            vlw_t entry(obt_tr.u64, obt_tr.u64 + n_words_per_obs);
            local_prob_histogram[entry]++;
        }

        shots_left -= shots_this_batch;
    }
    // Aggregate all the histogram results if using MPI.
    histogram_t prob_histogram;
    if (G_USE_MPI) {
        for (uint64_t r = 0; r < world_size; r++) {
            // Convert the histogram data into a vector (if rank == r).
            std::vector<std::pair<vlw_t, uint64_t>> hdata;
            if (world_rank == r) {
                for (auto pair : local_prob_histogram) {
                    hdata.push_back(pair);
                }
            }
            // Now, we will broadcast all the data in the histogram.
            //
            // First, broadcast the amount of data.
            uint n_histogram_entries = hdata.size();
            MPI_Bcast(&n_histogram_entries, 1, MPI_UNSIGNED, r, MPI_COMM_WORLD);
            for (uint i = 0; i < n_histogram_entries; i++) {
                vlw_t w;
                uint64_t cnt;
                if (world_rank == r) {
                    w = hdata[i].first;
                    cnt = hdata[i].second;
                }
                // Now, first broadcast the number of words in w.
                uint w_width = w.size();
                MPI_Bcast(&w_width, 1, MPI_UNSIGNED, r, MPI_COMM_WORLD);
                // Now, broadcast w itself.
                MPI_Bcast(&w[0], w_width, MPI_UNSIGNED_LONG, r, MPI_COMM_WORLD);
                // Finally, broadcast the count.
                MPI_Bcast(&cnt, 1, MPI_UNSIGNED_LONG, r, MPI_COMM_WORLD);
                // Update the final histogram.
                prob_histogram[w] += cnt;
            }
        }
    } else {
        prob_histogram = local_prob_histogram;
    }

    result_t res = {
        prob_histogram
    };
    return res;
}

void
SimManager::simulate_batch(uint64_t shots) {
    // Clear structures.
    sim->reset_sim();
    
    timing_table.clear();

    event_history_table.clear();
    obs_buffer_table.clear();
    pauli_frame_table.clear();

    meas_record_ctr = 0;
    max_obs_written = 0;

    // Execute batch.
    std::map<uint64_t, stim::simd_bits> pc_to_trials;
    // For ease of use, we will say that a shot is at a given
    // pc if the corresponding bit is 0.
    pc_to_trials[0] = stim::simd_bits(shots);
    pc_to_trials[0].clear();

    // br_data_t = (from_pc, to_pc, br_taken_pred)
    typedef std::tuple<uint64_t, uint64_t, stim::simd_bits> br_data_t;

    while (pc_to_trials.size()) {
        std::map<uint64_t, stim::simd_bits> next_pc_to_trials;
        std::vector<br_data_t> br_buffer;
        for (auto pair : pc_to_trials) {
            uint64_t pc = pair.first;
            stim::simd_bits shot_mask(pair.second);
            if (pc >= program.size())   continue;

            Instruction inst = program[pc];
            const std::vector<uint> operands(inst.operands);
            // Here, we will have a collection of variables that
            // are conditionally set based on the instruction.
            bool is_quantum_inst = ONLY_HAS_QUBIT_OPERANDS.count(inst.name);

            // Execute instructions.
            if (is_quantum_inst)    sim->snapshot();
            if (inst.name == "brifone" || inst.name == "brifzero") {
                const uint64_t br_pc = operands[0];
                const uint event = operands[1];
                stim::simd_bits br_taken(event_history_table[event]);
                if (inst.name == "brifzero")    update_mask.invert_bits();
                // Now we must filter out trials in br_taken that are
                // not at the current pc.
                br_taken.for_each_word(shot_mask,
                        [] (auto& brw, auto& smw)
                        {
                            // Keep br_taken = 1 wherever
                            // shot_mask = 0.
                            brw &= ~smw;
                            // Set shot_mask = 1 wherever
                            // br_taken = 1.
                            smw |= brw;
                        });
                br_buffer.push_back(std::tuple(pc, br_pc, br_taken));
            } else if (inst.name == "h") {
                sim->H(operands);
            } else if (inst.name == "x") {
                sim->X(operands);
            } else if (inst.name == "y") {
                sim->Y(operands);
            } else if (inst.name == "z") {
                sim->Z(operands);
            } else if (inst.name == "s") {
                sim->S(operands);
            } else if (inst.name == "measure") {
                sim->M(operands, params.errors.meas_record_ctr);
                meas_record_ctr += operands.size();
            } else if (inst.name == "reset") {
                sim->R(operands);
            } else if (inst.name == "cx") {
                sim->CX(operands);
            } else if (inst.name == "decode") {
                const uint frame = operands[0];
                auto eht_tr = event_history_table.transposed();
                for (uint64_t t = 0; t < shots; t++) {
                    if (G_FILTER_OUT_SYNDROMES
                        && eht_tr[t].popcnt() <= G_FILTERING_HAMMING_WEIGHT)
                    {
                        continue;
                    }
                    syndrome_t syndrome(eht_tr[t]);
                    auto res = decoder->decode_error(syndrome);
                    const uint obs = decoder->get_circuit().count_observables();
                    for (uint i = 0; i < obs; i++) {
                        pauli_frame_table[frame+i][t] ^= res.corr[i];
                    }
                }
            } else if (inst.name == "event") {
                const uint e = operands[0];
                for (uint i = 1; i < operands.size(); i++) {
                    uint x = operands[i];
                    event_history_table[e] ^= sim->record_table[x];
                }
            } else if (inst.name == "obs") {
                const uint obs = operands[0];
                for (uint i = 1; i < operands.size(); i++) {
                    uint x = operands[i];
                    obs_buffer_table[obs] ^= sim->record_table[x];
                }
                max_obs_written = max_obs_written < obs
                                    ? obs : max_obs_written;
            } else if (inst.name == "xorfr") {
                const uint obs = operands[0];
                const uint frame = operands[1];
                obs_buffer_table[obs] ^= pauli_frame_table[frame];
            }
            if (is_quantum_inst) {
                inject_operation_error(inst);
                timing_table[pc] += get_operation_latency(inst);
                if (!inst.annotations.count(Annotation::inject_timing_error)) {
                    inject_timing_error(timing_table[pc]);
                }

                sim->rollback_where(shot_mask);
            }
            next_pc_to_trials[pc+1] = shot_mask;
        }
        pc_to_trials = next_pc_to_trials;
    }
}

fp_t
SimManager::get_operation_latency(Instruction inst) {
    if (!IS_QUANTUM_INSTRUCTION.count(inst.name))   continue;

    if (inst.annotations.count(Annotation::no_tick))    return 0.0;

    std::vector<uint> operands(inst.get_qubit_operands());
    fp_t max_time_taken = 0.0;
    if (IS_2Q_OPERATION.count(inst.name)) {
        for (uint i = 0; i < operands.size(); i += 2) {
            uint x = operands[i];
            uint y = operands[i+1];
            auto x_y = std::make_pair(x, y);
            fp_t t = params.timing.op2q[inst.name][x_y];
            max_time_taken = t > max_time_taken ? t : max_time_taken;
        }
    } else {
        for (uint x : operands) {
            fp_t t = params.timing.op1q[inst.name][x];
            max_time_taken = t > max_time_taken ? t : max_time_taken;
        }
    }
    return max_time_taken;
}

void
SimManager::inject_operation_error(Instruction inst) {
    if (!IS_QUANTUM_INSTRUCTION.count(inst.name))   continue;

    bool inject_operation_error = !inst.annotations.count(Annotation::no_error);

    std::vector<uint> operands(inst.get_qubit_operands());
    if (IS_2Q_OPERATOR.count(inst.name)) {
        std::vector<fp_t> e_dp, e_li, e_lt;
        for (uint i = 0; i < operands.size(); i += 2) {
            uint x = operands[i];
            uint y = operands[i+1];
            auto x_y = std::make_pair(x, y);
            // Get error rates for (1) depolarizing errors,
            // (2) leakage errors, and (3) leakage transport errors.
            e_dp.push_back(params.errors.op2q[inst.name][x_y]);
            e_li.push_back(params.errors.op2q_leakage_injection[inst.name][x_y]);
            e_lt.push_back(params.errors.op2q_leakage_transport[inst.name][x_y]);
        }
        if (inject_operation_error) {
            sim->error_channel<StateSimulator::eLT>(operands, e_lt);
            sim->error_channel<StateSimulator::eLI>(operands, e_li);
            sim->error_channel<StateSimulator::eDP2>(operands, e_dp);
        }
    } else {
        fp_t e;
        for (uint x : operands) {
            e.push_back(params.errors.op1q[inst.name][x]);
        }
        if (inject_operation_error) {
            if (inst.name == "reset") {
                sim->error_channel<StateSimulator::eX>(operands, e);
            } else {
                sim->error_channel<StateSimulator::eDP1>(operands, e);
            }
        }
    }
}

void
SimManager::inject_timing_error(fp_t time) {
    for (uint i = 0; i < n_qubits; i++) {
        fp_t t1 = params.timing.t1[i];
        fp_t t2 = params.timing.t2[i];
        fp_t e_ad = 0.25*(1 - exp(-time/t1));
        fp_t e_pd = 0.5*(1 - exp(-time/t2));
        // Amplitude damping is a non-Clifford channel.
        //
        // Current implementation (here) is via Pauli twirling:
        sim->error_channel<StateSimulator::eX>(e_ad);
        sim->error_channel<StateSimulator::eY>(e_ad);
        sim->error_channel<StateSImulator::eZ>(e_pd - e_ad);
    }
}

}   // qontra
