/* author: Suhas Vittal
 *  date:   10 December 2023
 * */

#include "defs/filesystem.h"
#include "experiments.h"
#include "instruction.h"
#include "sim/memory_sim.h"

#include <fstream>
#include <iostream>

#include <math.h>
#include <mpi.h>
#include <stdio.h>

namespace qontra {

using namespace graph;
using namespace lattice;

MemorySimulator::MemorySimulator(LatticeGraph& gr)
    :// Global variables
    s_lrcs_per_round(MPI_SUM, 1),
    s_leakage_population_ratio(MPI_SUM, 32),
    // Private variables
    lattice_graph(gr),
    meas_ctr(0),
    event_ctr(0),
    obs_ctr(0),
    elapsed_time(0),
    shot_time_delta_map(),
    sample_circuit(),
    is_recording_stim_instructions(false),
    // Dummy initializations:
    all_qubits(),
    data_qubits(),
    parity_qubits(),
    zp_qubits(),
    xp_qubits(),
    syndrome_parity_qubits(),
    obs_list(),
    n_qubits(0),
    n_detection_events(0),
    n_observables(0),
    sim(nullptr),
    syndromes(1, 1),
    // Variables for LRC:
    lrc_await_queue(),
    lrc_optimal_lrc_map_table(),
    // Variables for ERASER:
    eraser_qubit_to_syndrome_buffer_index_map(),
    eraser_swap_lookup_table(),
    eraser_recently_scheduled_qubits_table(1,1),
    eraser_syndrome_buffer(1,1)
{
    std::vector<sptr<vertex_t>> vertex_list = lattice_graph.get_vertices();
    for (sptr<vertex_t> v : vertex_list) {
        const uint64_t id = v->id;
        all_qubits.push_back(id);
        if (v->qubit_type == vertex_t::type::data) {
            data_qubits.push_back(id);
        } else {
            parity_qubits.push_back(id);
            if (v->qubit_type == vertex_t::type::zparity) {
                zp_qubits.push_back(id);
            } else {
                xp_qubits.push_back(id);
            }
        }
    }
    n_qubits = all_qubits.size();
    sim = std::make_unique<FrameSimulator>(n_qubits, experiments::G_SHOTS_PER_BATCH);
    reset();
}

void
MemorySimulator::initialize() {
    // Clear out the stats.
    s_leakage_population_ratio.resize(config.rounds);

    s_lrcs_per_round.fill(0);
    s_leakage_population_ratio.fill(0);

    // We'll just remake the following data in case the config has changed.
    syndrome_parity_qubits = config.is_memory_x ? xp_qubits : zp_qubits;
    auto vertex_obs_list = config.is_memory_x ? lattice_graph.x_obs_list : lattice_graph.z_obs_list;
    obs_list.clear();
    for (std::vector<sptr<vertex_t>> arr : vertex_obs_list) {
        std::vector<uint> obs;
        for (sptr<vertex_t> v : arr)  obs.push_back(v->id);
        obs_list.push_back(obs);
    }

    n_detection_events = (config.rounds+1) * syndrome_parity_qubits.size();
    n_observables = obs_list.size();

    syndromes = stim::simd_bit_table(n_detection_events+n_observables, experiments::G_SHOTS_PER_BATCH);

#ifdef QONTRA_MEMORY_SIM_EXT_ENABLED
    eraser_initialize();
#endif
}

void
MemorySimulator::reset() {
    meas_ctr = 0;
    event_ctr = 0;
    obs_ctr = 0;
    elapsed_time = 0;
    shot_time_delta_map.clear();

    sim->reset_sim();
    syndromes.clear();

#ifdef QONTRA_MEMORY_SIM_EXT_ENABLED
    lrc_reset();
    eraser_reset();
#endif
}

MemorySimulator::time_t
MemorySimulator::do_gate(std::string op, std::vector<uint> operands, int64_t trial) {
    // Note -- we do not want to record trial-specific errors and instructions.
    // Thus, we only record instructions/errors in Stim if trial < 0.
    if (op == "h") {
        sim->H(operands, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("H", operands);
        }
    } else if (op == "x") {
        sim->X(operands, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("X", operands);
        }
    } else if (op == "z") {
        sim->Z(operands, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("Z", operands);
        }
    } else if (op == "cx") {
        sim->CX(operands, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("CX", operands);
        }
    } else if (op == "reset") {
        sim->R(operands, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("R", operands);
        }
    } else if (op == "liswap") {
        sim->LEAKAGE_ISWAP(operands, trial);
    }
    else return 0;  // Nothing to be done -- treat as NOP.
    // Inject any errors.
    bool is_2q = IS_2Q_OPERATOR.count(op); 
    time_t instruction_latency = 0.0;   // The latency of the SIMD operation is the maximum
                                        // latency of any sub-operation.
    if (is_2q) {
        std::vector<fp_t> dp_array;
        std::vector<fp_t> li_array;
        std::vector<fp_t> lt_array;
        for (uint i = 0; i < operands.size(); i += 2) {
            uint q1 = operands[i], q2 = operands[i+1];
            auto q1_q2 = std::make_pair(q1, q2);

            fp_t e_dp = config.errors.op2q[op][q1_q2];
            fp_t e_li = config.errors.op2q_leakage_injection[op][q1_q2];
            fp_t e_lt = config.errors.op2q_leakage_transport[op][q1_q2];
            
            lt_array.push_back(e_lt);
            li_array.push_back(e_li);
            dp_array.push_back(e_dp);
            if (is_recording_stim_instructions && trial < 0) {
                if (e_dp > 0) sample_circuit.safe_append_ua("DEPOLARIZE2", {q1, q2}, e_dp);
            }
            // Get instruction latency.
            fp_t t = config.timing.op2q[op][q1_q2];
            if (t > instruction_latency) instruction_latency = t;
        }
        if (trial >= 0) {
            // We are injecting errors on a specific trial.
            for (uint i = 0; i < lt_array.size(); i++) {
                uint q1 = operands[2*i], q2 = operands[2*i+1];
                if (sim->get_probability_sample_from_rng() < lt_array[i]) sim->eLT(q1, q2, trial);
                if (sim->get_probability_sample_from_rng() < li_array[i]) sim->eLI(q1, q2, trial);
                if (sim->get_probability_sample_from_rng() < dp_array[i]) sim->eDP2(q1, q2, trial);
            }
        } else {
            sim->error_channel<&StateSimulator::eLT>(operands, lt_array);
            sim->error_channel<&StateSimulator::eLI>(operands, li_array);
            sim->error_channel<&StateSimulator::eDP2>(operands, dp_array);
        }
    } else {
        std::vector<fp_t> error_rates;
        for (uint i = 0; i < operands.size(); i++) {
            uint q = operands[i];
            fp_t e = config.errors.op1q[op][q];

            error_rates.push_back(e);

            if (is_recording_stim_instructions && trial < 0) {
                if (e > 0) {
                    if (op == "reset") {
                        sample_circuit.safe_append_ua("X_ERROR", {q}, e);
                    } else {
                        sample_circuit.safe_append_ua("DEPOLARIZE1", {q}, e);
                    }
                }
            }
            // Get instruction latency.
            fp_t t = config.timing.op1q[op][q];
            if (t > instruction_latency) instruction_latency = t;
        }
        if (trial >= 0) {
            // We are injecting errors for a specific trial.
            for (uint i = 0; i < error_rates.size(); i++) {
                uint q = operands[i];
                if (sim->get_probability_sample_from_rng() < error_rates[i]) {
                    if (op == "reset")  sim->eX(q, trial);
                    else                sim->eDP1(q, trial);
                }
            }
        } else {
            if (op == "reset") {
                sim->error_channel<&StateSimulator::eX>(operands, error_rates);
            } else {
                sim->error_channel<&StateSimulator::eDP1>(operands, error_rates);
            }
        }
    }
    return instruction_latency;
}

MemorySimulator::time_t
MemorySimulator::do_measurement(std::vector<uint> operands, int64_t trial) {
    // Like with do_gate, we will only record instructions for Stim if trial < 0,
    // as this implies a measurement done across all trials.
    //
    // TODO: implement non-deterministic measurement counters.
    //
    std::vector<fp_t> m1w0_array, m0w1_array;
    time_t instruction_latency = 0.0;
    for (uint i = 0; i < operands.size(); i++) {
        uint q = operands[i];
        m1w0_array.push_back(config.errors.m1w0[q]);
        m0w1_array.push_back(config.errors.m0w1[q]);
        // Stim doesn't have a way to implement biased measurements.
        // So, we just take the mean.
        if (is_recording_stim_instructions && trial < 0) {
            fp_t emean = 0.5 * (config.errors.m1w0[q]+config.errors.m0w1[q]);
            sample_circuit.safe_append_ua("X_ERROR", {q}, emean);
        }
        // Update measurement counter (only if trial < 0 -- otherwise we run into issues).
        if (trial < 0)  meas_ctr_map[q] = meas_ctr+i;
        // Get instruction latency.
        fp_t t = config.timing.op1q["measure"][q];
        if (t > instruction_latency) instruction_latency = t;
    }
    // Perform the measurement.
    sim->M(operands, m1w0_array, m0w1_array, meas_ctr, trial);
    if (is_recording_stim_instructions && trial < 0) {
        sample_circuit.safe_append_u("M", operands);
    }
    // We will only update the counter if this is a common measurement.
    if (trial < 0)  meas_ctr += operands.size();

    return instruction_latency;
}

void
MemorySimulator::inject_idling_error_positive(std::vector<uint> on_qubits, int64_t trial) {
    // Do NOT record the error if trial >= 0.
    std::vector<fp_t> error_rates;
    for (uint q : on_qubits) {
        fp_t e = config.errors.idling[q];
        error_rates.push_back(e);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_ua("DEPOLARIZE1", {q}, e);
        }
    }

    if (trial >= 0) {
        // We are only injecting the error on a single trial.
        for (size_t i = 0; i < error_rates.size(); i++) {
            uint q = on_qubits[i];
            if (sim->get_probability_sample_from_rng() < error_rates[i]) sim->eDP1(q, trial);
        }
    } else {
        sim->error_channel<&StateSimulator::eDP1>(on_qubits, error_rates);
    }
}

void
MemorySimulator::inject_idling_error_negative(std::vector<uint> not_on_qubits, int64_t trial) {
    std::vector<uint> on_qubits;
    for (size_t i = 0; i < n_qubits; i++) {
        if (std::find(not_on_qubits.begin(), not_on_qubits.end(), i) == not_on_qubits.end()) {
            on_qubits.push_back(i);
        }
    }
    inject_idling_error_positive(on_qubits, trial);
}

void
MemorySimulator::create_event_or_obs(std::vector<uint> operands, bool create_event) {
    uint64_t index = create_event ? event_ctr : n_detection_events+obs_ctr;
    std::vector<uint> offsets;   // Needed for Stim circuit.
    for (uint meas_time : operands) {
        syndromes[index] ^= sim->record_table[meas_time];
        int diff = meas_ctr - meas_time;
        offsets.push_back(stim::TARGET_RECORD_BIT | diff);
    }

    if (is_recording_stim_instructions) {
        if (create_event) {
            sample_circuit.safe_append_u("DETECTOR", offsets);
        } else {
            sample_circuit.safe_append_u("OBSERVABLE_INCLUDE", offsets, {static_cast<fp_t>(obs_ctr)});
        }
    }
    if (create_event)   event_ctr++;
    else                obs_ctr++;
}

void
MemorySimulator::inject_timing_error(std::vector<uint> qubits) {
    std::vector<fp_t> xy_array, z_array, li_array;
    for (uint q : qubits) {
        fp_t t1 = config.timing.t1[q];
        fp_t t2 = config.timing.t2[q];
        
        fp_t e_ad = 0.25 * (1 - exp(-elapsed_time/t1));
        fp_t e_pd = 0.5 * (1 - exp(-elapsed_time/t2));

        xy_array.push_back(e_ad);
        z_array.push_back(e_pd-e_ad);
        li_array.push_back(0.1*e_ad);

        if (is_recording_stim_instructions) {
            sample_circuit.safe_append_ua("X_ERROR", {q}, e_ad);
            sample_circuit.safe_append_ua("Y_ERROR", {q}, e_ad);
            sample_circuit.safe_append_ua("Z_ERROR", {q}, e_pd-e_ad);
        }
    }
    if (config.enable_leakage) {
        sim->error_channel<&StateSimulator::eL>(qubits, li_array);
    }
    sim->error_channel<&StateSimulator::eX>(qubits, xy_array);
    sim->error_channel<&StateSimulator::eY>(qubits, xy_array);
    sim->error_channel<&StateSimulator::eZ>(qubits, z_array);
    // If there are any trials with time deltas, handle them now.
    for (auto pair : shot_time_delta_map) {
        uint64_t shot = pair.first;
        time_t delta = pair.second;
        for (uint q : qubits) {
            fp_t t1 = config.timing.t1[q];
            fp_t t2 = config.timing.t2[q];
            
            fp_t e_ad = 0.25 * (1 - exp(-delta/t1));
            fp_t e_pd = 0.5 * (1 - exp(-delta/t2));

            if (sim->get_probability_sample_from_rng() < 0.1*e_ad)   sim->eL(q, shot);
            if (sim->get_probability_sample_from_rng() < e_ad)       sim->eX(q, shot);
            if (sim->get_probability_sample_from_rng() < e_ad)       sim->eY(q, shot);
            if (sim->get_probability_sample_from_rng() < e_pd-e_ad)  sim->eZ(q, shot);
        }
    }
    elapsed_time = 0;
    shot_time_delta_map.clear();
}

void
MemorySimulator::run(uint64_t shots) {
    using namespace experiments;

    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }
    const uint64_t local_shots = shots / world_size + (world_rank == 0 ? shots % world_size : 0);

    uint64_t shots_left = local_shots;
    int batchno = world_rank;

    // Set the seed. For the first batch only, the 0th process will record instructions to
    // write to a Stim file.
    is_recording_stim_instructions = (world_rank == 0);
    sim->set_seed(G_BASE_SEED + world_rank);

    initialize();

    // Synchronize before starting.
    if (G_USE_MPI)  MPI_Barrier(MPI_COMM_WORLD);

    while (shots_left) {
        const uint64_t shots_this_batch = shots_left < G_SHOTS_PER_BATCH ? shots_left : G_SHOTS_PER_BATCH;
        run_batch(shots_this_batch);
        // Record results.
        if (is_recording_stim_instructions && world_rank == 0) {
            std::ofstream stim_out(config.stim_output_file);
            stim_out << sample_circuit.str() << "\n";
            is_recording_stim_instructions = false;
        }
        // Record syndromes.
        std::string syndrome_output_file = config.syndrome_output_folder 
                                            + "/batch_"
                                            + std::to_string(batchno)
                                            + ".dets";
        FILE* fout = fopen(syndrome_output_file.c_str(), "w");
        const uint64_t n_syndrome_bits = n_detection_events + n_observables;
        const stim::simd_bits ref(n_syndrome_bits);
        stim::write_table_data(fout, shots_this_batch, n_syndrome_bits, ref, syndromes,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS, 'D', 'L', n_detection_events);
        fclose(fout);

        batchno += world_size;
        shots_left -= shots_this_batch;
    }
    // 
    // Accumulate statistics.
    //
    s_lrcs_per_round.reduce();
    s_leakage_population_ratio.reduce();

    s_lrcs_per_round /= static_cast<fp_t>(config.rounds * shots);
    s_leakage_population_ratio /= static_cast<fp_t>(n_qubits * shots);
    // Write statistics to output file.
    if (world_rank == 0) {
        bool write_header = file_exists(config.data_output_file);
        safe_create_directory(config.data_output_file.c_str());

        std::ofstream fout(config.data_output_file, std::ios::app);
        if (write_header) {
            // The file does not exist.

            fout << "Distance,"
                    << "Rounds,"
                    << "CX Error Rate,"
                    << "LRC Policy,"
                    << "LRC Circuit,"
                    << "LRCs per Round";
            for (uint r = 0; r < config.rounds; r++) {
                fout << ",LPR Round " << r;
            }
            fout << "\n";
        }
        // Write stats.
        std::string lrc_policy_name;
        if (config.lrc_policy == lrc_policy_t::none)            lrc_policy_name = "N/A";
        else if (config.lrc_policy == lrc_policy_t::always)     lrc_policy_name = "Always-LRCs";
        else if (config.lrc_policy == lrc_policy_t::optimal)    lrc_policy_name = "Optimal-LRCs";
        else if (config.lrc_policy == lrc_policy_t::eraser)     lrc_policy_name = "ERASER";

        std::string lrc_circuit_name;
        if (config.lrc_circuit == lrc_circuit_t::swap)  lrc_circuit_name = "SWAP";
        if (config.lrc_circuit == lrc_circuit_t::dqlr)  lrc_circuit_name = "DQLR";

        fout << config.distance << ","
                << config.rounds << ","
                << config.errors.op2q["cx"][std::make_pair(0, 1)] << ","
                << lrc_policy_name << ","
                << lrc_circuit_name << ","
                << s_lrcs_per_round.at();
        for (uint r = 0; r < config.rounds; r++) {
            fout << "," << s_leakage_population_ratio.at(r);
        }
        fout << "\n";
    }
}

void
MemorySimulator::run_batch(uint64_t shots) {
    sim->shots = shots;
    const bool mx = config.is_memory_x;

    reset();
    //
    // PROLOGUE
    //
    elapsed_time += do_gate("reset", all_qubits);
    if (mx) {
        elapsed_time += do_gate("h", data_qubits);
    }
    //
    // BODY
    //
    for (uint r = 0; r < config.rounds; r++) {
        std::map<uint, uint64_t> prev_meas_ctr_map(meas_ctr_map);

#ifdef QONTRA_MEMORY_SIM_EXT_ENABLED
        stim::simd_bits lrc_shots_with_leakage = lrc_optimal_identify_lrcs();
#endif

        inject_timing_error(data_qubits);
        elapsed_time += do_gate("h", xp_qubits);
        inject_idling_error_negative(xp_qubits);

        uint cx_depth = 0;
        while (1) {
            std::vector<uint> cx_operands;
            for (uint pi : parity_qubits) {
                auto pv = lattice_graph.get_vertex(pi);
                for (auto dv : lattice_graph.get_neighbors(pv)) {
                    auto e = lattice_graph.get_edge(pv, dv);
                    if (e->cx_time != cx_depth) continue;
                    uint di = dv->id;
                    if (pv->qubit_type == vertex_t::type::xparity) {
                        cx_operands.push_back(pi);
                        cx_operands.push_back(di);
                    } else {
                        cx_operands.push_back(di);
                        cx_operands.push_back(pi);
                    }
                    break;
                }
            }
            if (cx_operands.empty()) break;
            elapsed_time += do_gate("cx", cx_operands);
            inject_idling_error_negative(cx_operands);
            cx_depth++;
        }
        elapsed_time += do_gate("h", xp_qubits);
        inject_idling_error_negative(xp_qubits);

#ifdef QONTRA_MEMORY_SIM_EXT_ENABLED
        stim::simd_bits eraser_shots_using_lrcs(1);
        if (config.lrc_policy == lrc_policy_t::always) {
            lrc_execute_lrcs_from_await_queue();
            goto memory_sim_make_detection_events;
        } else if (config.lrc_policy == lrc_policy_t::optimal) {
            lrc_optimal_perform_lrcs();
            // Save the state of the simulator.
            sim->snapshot();
        } else if (config.lrc_policy == lrc_policy_t::eraser) {
            eraser_shots_using_lrcs = eraser_execute_lrcs();
            sim->snapshot();
        }
#endif

        elapsed_time += do_measurement(parity_qubits);
        inject_idling_error_positive(data_qubits);
        elapsed_time += do_gate("reset", parity_qubits);

#ifdef QONTRA_MEMORY_SIM_EXT_ENABLED
        if (config.lrc_policy == lrc_policy_t::optimal) {
            // Now, rollback any changes on trials that had LRCs.
            sim->rollback_where(lrc_shots_with_leakage);
        } else if (config.lrc_policy == lrc_policy_t::eraser) {
            sim->rollback_where(eraser_shots_using_lrcs);
            // Update Eraser's syndrome buffer.
            eraser_update_syndrome_buffer(prev_meas_ctr_map);
        }
#endif

memory_sim_make_detection_events:
        // Create detection events.
        for (uint pi : syndrome_parity_qubits) {
            std::vector<uint> meas_list{ (uint) meas_ctr_map[pi] };
            if (r > 0) {
                meas_list.push_back(prev_meas_ctr_map[pi]);
            }
            create_event_or_obs(meas_list);
        }
        //
        // Update statistics:
        //      s_leakage_population_ratio
        //
        for (uint q = 0; q < n_qubits; q++) {
            s_leakage_population_ratio[r] += static_cast<fp_t>(sim->leak_table[q].popcnt());
        }
    }
    //
    // EPILOGUE
    //
    if (mx) do_gate("h", data_qubits);
    do_measurement(data_qubits);
    for (uint pi : syndrome_parity_qubits) {
        std::vector<uint> meas_list{ (uint) meas_ctr_map[pi] };
        auto pv = lattice_graph.get_vertex(pi);
        for (auto dv : lattice_graph.get_neighbors(pv)) {
            uint di = dv->id;
            meas_list.push_back(meas_ctr_map[di]);
        }
        create_event_or_obs(meas_list);
    }
    for (auto obs_qubits : obs_list) {
        std::vector<uint> meas_list;
        for (uint di : obs_qubits) meas_list.push_back(meas_ctr_map[di]);
        create_event_or_obs(meas_list, false);
    }
}

}   // qontra
