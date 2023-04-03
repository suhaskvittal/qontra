/*
 *  author: Suhas Vittal
 *  date:   9 August 2022
 * */

#include "benchmark.h"

namespace qrc {

#define STATBENCH_DEBUG
//#define STATBENCH_DEBUG2

/* Helper Functions */
inline fp_t __CHS(fp_t x, fp_t y, fp_t z) {
    return z >= 0 ? z : (y >=0 ? y : x);
}

/* Main Functions */

fp_t
min(const std::vector<fp_t>& data) {
    return *std::min_element(data.begin(), data.end());
}

fp_t
max(const std::vector<fp_t>& data) {
    return *std::max_element(data.begin(), data.end());
}

fp_t
mean(const std::vector<fp_t>& data) {
    fp_t sum = std::accumulate(data.begin(), data.end(), 0);
    return sum / data.size();
}
fp_t
stdev(const std::vector<fp_t>& data) {
    fp_t sum = 0.0;
    fp_t m = mean(data);
    for (uint32_t i = 0; i < data.size(); i++) {
        sum += (data[i] - m)*(data[i] - m);
    }
    return sqrt(sum/data.size());
}

void
b_decoder_ler(Decoder * decoder_p, uint64_t shots, std::mt19937_64& rng,
        bool save_per_shot_data) 
{
    const uint64_t _shots = shots;
    // Clear stats.
    decoder_p->clear_stats();
    // Declare statistics
    uint32_t array_size = save_per_shot_data ? shots : 1;

    uint64_t total_shots = shots;
    std::vector<std::vector<uint8_t>> syndromes(array_size);
    std::vector<fp_t> execution_times(array_size);
    std::vector<fp_t> memory_overheads(array_size);
    uint32_t n_logical_errors = 0;
    fp_t mean_execution_time = 0.0,
         max_execution_time = 0.0,
         max_execution_time_for_correctable = 0.0;

    std::map<uint, uint64_t> hamming_weight_freq;

    uint32_t sn = 0;
    uint32_t bn = 0;
    while (shots > 0) {
        uint32_t shots_this_round = shots > MAX_SHOTS ? MAX_SHOTS : shots;
        // Get samples from Stim.
        stim::simd_bit_table leakage_buffer(1, 1);
        stim::simd_bit_table sample_buffer = 
            stim::detector_samples(decoder_p->circuit, shots_this_round,
                    false, true, rng, true, leakage_buffer);
        // Last part of samples is the actual observable.
        // We are trying to match that.
        sample_buffer = sample_buffer.transposed();
        leakage_buffer = leakage_buffer.transposed();
        uint n_detectors = decoder_p->circuit.count_detectors();
        uint n_observables = decoder_p->circuit.count_observables();

        for (uint32_t i = 0; i < shots_this_round; i++) {
            auto syndrome = 
                _to_vector(sample_buffer[i], n_detectors, n_observables);
            uint hw = 
                std::accumulate(syndrome.begin(), syndrome.end() - n_observables, 0);
            if (hw & 0x1) {
                hw++;
            }
            if (!hamming_weight_freq.count(hw)) {
                hamming_weight_freq[hw] = 0;
            }
            hamming_weight_freq[hw]++;
            // Update some stats before decoding.
            if (save_per_shot_data) {
                syndromes[sn] = syndrome;
            }
            if (hw > 0) {
                DecoderShotResult res = decoder_p->decode_error(syndrome);
                // Update statistics.
//                n_logical_errors += (res.is_logical_error || leakage_buffer[i][n_detectors]);
                n_logical_errors += res.is_logical_error;
                mean_execution_time += res.execution_time / ((fp_t)total_shots);
                if (res.execution_time > max_execution_time) {
                    max_execution_time = res.execution_time;
                    if (!res.is_logical_error) {
                        max_execution_time_for_correctable = res.execution_time;    
                    }
                }
                if (save_per_shot_data) {
                    execution_times[sn] = res.execution_time;
                    memory_overheads[sn] = res.memory_overhead;
                }
            } else {
                if (save_per_shot_data) {
                    execution_times[sn] = 0;
                    memory_overheads[sn] = 0;
                }
            }
            sn++;
        }
        shots -= shots_this_round;
        bn++;
    }
    // Update stats in decoder.
    decoder_p->syndromes = syndromes;
    decoder_p->execution_times = execution_times;
    decoder_p->memory_overheads = memory_overheads;
    decoder_p->n_logical_errors = n_logical_errors;
    decoder_p->mean_execution_time = mean_execution_time;
    decoder_p->max_execution_time = max_execution_time;
    decoder_p->max_execution_time_for_correctable = 
        max_execution_time_for_correctable;
    decoder_p->hamming_weight_dist.clear();
    for (auto pair : hamming_weight_freq) {
        decoder_p->hamming_weight_dist[pair.first] = ((fp_t) pair.second) / ((fp_t) _shots);
    }
}

benchmark::StatisticalResult
b_statistical_ler(
    dgf_t& mkdec, 
    uint code_dist,
    fp_t p_ref,
    uint64_t shots_per_call, 
    std::mt19937_64& rng, 
    bool use_mpi)
{
    int world_rank = 0;
    int world_size = 1;
    if (use_mpi) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }

    Decoder * curr_decoder = mkdec(p_ref);

    benchmark::StatisticalResult statres;
    fp_t p_curr = p_ref;
    while (statres.n_logical_errors < 100) {
        uint64_t shots = shots_per_call;
        uint64_t local_shots = shots / world_size;
        if (world_rank == 0) {
            local_shots += shots % world_size;
        }

        uint64_t local_errors = 0;

        stim::FrameSimulator sim(circuit.count_qubits(), shots_this_round, SIZE_MAX, rng);
        sim.reference_error_rate = p_ref;
        while (local_shots) {
            uint64_t shots_this_round = local_shots < MAX_SHOTS ? local_shots : MAX_SHOTS;
            const stim::Circuit circ(curr_decoder->circuit);
            sim.reset_all_and_run(circ);
            
            stim::simd_bit_table result_table(1, 1);
            stim::simd_bit_table leakage_table(1, 1);  // Unused.
            stim::read_from_sim(
                    sim,
                    stim::DetectorsAndObservables(circ),
                    false,
                    true,
                    false,
                    result_table,
                    leakage_table);
            result_table = result_table.transposed();
            for (uint64_t s = 0; s < shots_this_round; s++) {
                auto syndrome = _to_vector(result_table[s],
                                            circ.count_detectors(),
                                            circ.count_observables());
                auto decoder_res = curr_decoder->decode_error(syndrome);
                if (decoder_res.is_logical_error) {
                    local_errors++;
                }
            }
            
            local_shots -= shots_this_round;
        }

        fp_t total_log_prob_sim, total_log_prob_ref;
        if (use_mpi) {
            MPI_Allreduce(&sim.log_prob_sim, &total_log_prob_sim, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&sim.log_prob_reference, &total_log_prob_ref, 1, MPI_DOUBLE,
                            MPI_SUM, MPI_COMM_WORLD);
        } else {
            total_log_prob_sim = sim.log_prob_sim;
            total_log_prob_ref = sim.log_prob_reference;
        }

        fp_t stat_shots = shots * pow(10, total_log_prob_sim - total_log_prob_ref);
        statres.true_shots += shots;
        statres.statistical_shots += stat_shots;

        uint64_t n_logical_errors;
        MPI_Allreduce(&local_errors, &n_logical_errors, 1, MPI_UNSIGNED_LONG,
                            MPI_SUM, MPI_COMM_WORLD);
        statres.n_logical_errors += n_logical_errors;

        p_curr *= 2;
        curr_decoder = mkdec(p_curr);
    }
    delete curr_decoder;

    return statres;
}

stim::Circuit
build_circuit(
    uint code_dist,
    fp_t error_mean,
    fp_t error_stddev,
    bool is_memory_z,
    bool is_rotated,
    bool both_stabilizers,
    uint8_t other_flags,
    uint rounds,
    fp_t clevel_error_mean,
    fp_t clevel_error_stddev,
    fp_t pauliplus_error_mean,
    fp_t pauliplus_error_stddev,
    fp_t round_dp_mean,
    fp_t sq_dp_mean,
    fp_t cx_dp_mean,
    fp_t reset_flip_mean,
    fp_t meas_flip_mean,
    fp_t round_dp_stddev,
    fp_t sq_dp_stddev,
    fp_t cx_dp_stddev,
    fp_t reset_flip_stddev,
    fp_t meas_flip_stddev,
    fp_t round_leak_mean,
    fp_t clifford_leak_mean,
    fp_t leak_transport_mean,
    fp_t round_leak_stddev,
    fp_t clifford_leak_stddev,
    fp_t leak_transport_stddev)
{
    if (rounds == 0) {
        rounds = code_dist;
    }
    std::string circ_type = (is_rotated ? "" : "un");
    circ_type += "rotated_memory_";
    circ_type += (is_memory_z ? "z" : "x");

    stim::CircuitGenParameters params(rounds, code_dist, circ_type);
    // Declare error rates.
    params.before_round_data_depolarization = __CHS(error_mean, clevel_error_mean, round_dp_mean);
    params.after_clifford_depolarization = __CHS(error_mean, clevel_error_mean, cx_dp_mean);
    params.after_clifford_sq_depolarization = __CHS(error_mean, clevel_error_mean, sq_dp_mean);
    params.after_reset_flip_probability = __CHS(error_mean, clevel_error_mean, reset_flip_mean);
    params.before_measure_flip_probability = __CHS(error_mean, clevel_error_mean, meas_flip_mean);

    params.before_round_data_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, round_dp_stddev);
    params.after_clifford_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, cx_dp_stddev);
    params.after_clifford_sq_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, sq_dp_stddev);
    params.after_reset_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, reset_flip_stddev);
    params.before_measure_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, meas_flip_stddev);

    params.before_round_leakage_probability = __CHS(error_mean, pauliplus_error_mean, round_leak_mean);
    params.after_clifford_leakage_probability = __CHS(error_mean, pauliplus_error_mean, clifford_leak_mean);
    params.after_clifford_leakage_transport = __CHS(error_mean, pauliplus_error_mean, leak_transport_mean);
    
    params.before_round_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, round_leak_stddev);
    params.after_clifford_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, clifford_leak_stddev);
    params.after_clifford_leakage_transport_stddev = __CHS(error_stddev, pauliplus_error_stddev, leak_transport_stddev);

    params.both_stabilizers = both_stabilizers;
    params.swap_lru = other_flags & BC_FLAG_SWAP_LRU_V1;
    params.swap_lru_with_no_swap = other_flags & BC_FLAG_SWAP_LRU_V2;
    params.initial_state_is_basis_1 = other_flags & BC_FLAG_INVERT_STATE;

    stim::Circuit circ = generate_surface_code_circuit(params).circuit;
    return circ;
}

}  // qrc
