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
            DecoderShotResult res = decoder_p->decode_error(syndrome);
            // Update statistics.
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
    fp_t start_p, 
    fp_t final_p,
    uint64_t shots, 
    uint64_t update_rate,
    std::mt19937_64& rng, 
    bool use_mpi,
    bool bootstrap_model,
    std::map<uint, uint64_t> bootstrap_data,
    fp_t use_bootstrap_model_until_p)
{
    int world_rank = 0;
    int world_size = 1;
    if (use_mpi) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }

    Decoder * curr_decoder = mkdec(start_p);
    benchmark::ErrorDistribution * error_dist;
    benchmark::ErrorDistribution * alt_dist;
    if (bootstrap_model) {
        auto numerical_dist = benchmark::statbench::build_numerical_distribution(
                                    bootstrap_data,
                                    benchmark::statbench::NumericalModel::normal);
        auto analytical_dist = benchmark::statbench::bootstrap_from_normal_approximation(
                                    curr_decoder->circuit, 
                                    benchmark::statbench::AnalyticalModel::poisson,
                                    numerical_dist);
        error_dist = &numerical_dist;
        alt_dist = &analytical_dist;
    } else {
        auto analytical_dist = benchmark::statbench::build_analytical_distribution(
                                    curr_decoder->circuit, benchmark::statbench::AnalyticalModel::binomial);
        error_dist = &analytical_dist;
        alt_dist = &analytical_dist;
    }

    const fp_t log_p_update = ((fp_t)update_rate) * (log(final_p) - log(start_p)) / ((fp_t) shots);
    const fp_t log_mean_prob = error_dist->pmf(2*error_dist->get_mean(), true);

#ifdef STATBENCH_DEBUG
    if (world_rank == 0) {
        std::cout << "expected hamming weight: " << error_dist->get_mean() 
                    << "\tvariance: " << error_dist->get_variance()
                    << "\tlog(prob) = " << log_mean_prob << "\n";
        std::cout << "log_p_update = " << log_p_update << "\n";
    }
#endif

    benchmark::StatisticalResult statres;
    statres.true_shots = shots;

    std::map<uint, fp_t> hamming_weight_freq;

    fp_t pc = start_p;    // Current physical error rate
    fp_t prev_sample_ss = shots < update_rate ? shots : update_rate;

    bool first_iter = true;  // We need to adjust any inaccuracies in our model in the first iteration only.
    while (shots > 0) {
#ifdef STATBENCH_DEBUG
        if (world_rank == 0) {
            std::cout << "shots left: " << shots << "\n";
        }
#endif
        const uint64_t shots_this_round = shots < update_rate ? shots : update_rate;
        // Benchmark decoder
        uint64_t shots_this_sample = shots_this_round / world_size;
        if (world_rank == world_size - 1) {
            shots_this_sample += shots_this_round % world_size;
        }

        std::map<uint, fp_t> sample_hw_freq;
        uint64_t sample_errors;
        fp_t sample_statistical_shots;
        fp_t sample_mean;

        std::map<uint, fp_t> local_freq;
        uint64_t local_errors = 0;
        fp_t local_statistical_shots = 0.0;
        fp_t hwsum = 0;
        while (shots_this_sample > 0) {
            uint64_t shots_this_batch = shots_this_sample < 100'000 ? shots_this_sample : 100'000;
            stim::simd_bit_table sample_buffer = 
                stim::detector_samples(curr_decoder->circuit, shots_this_batch,
                        false, true, rng);
            sample_buffer = sample_buffer.transposed();
            for (uint64_t i = 0; i < shots_this_batch; i++) {
                std::vector<uint8_t> syndrome = _to_vector(sample_buffer[i], 
                                                        curr_decoder->circuit.count_detectors(),
                                                        curr_decoder->circuit.count_observables());
                uint hw = std::accumulate(syndrome.begin(),
                                        syndrome.begin() + curr_decoder->circuit.count_detectors(),
                                        0);
                if (hw & 0x1) {
                    hw++;
                }
                hwsum += hw;
                if (!local_freq.count(hw)) {
                    local_freq[hw] = 0;
                }
                local_freq[hw]++;
                if (hw) {
                    auto res = curr_decoder->decode_error(syndrome);
                    local_errors += res.is_logical_error;
                }
                // Compute statistical shots for this batch.
                fp_t log_ss;
                fp_t log_hw_prob = error_dist->pmf(hw, true);
                if (-log_hw_prob > 2.5+log(statres.statistical_shots)) {
                    log_hw_prob = alt_dist->pmf(hw, true);
                }
                log_ss = log_mean_prob - log_hw_prob;
                local_statistical_shots += pow(M_E, log_ss);
            }
            
            shots_this_sample -= shots_this_batch;
        }
        // Collect results.
        for (uint i = 0; i < curr_decoder->circuit.count_detectors(); i++) {
            fp_t f = 0.0;
            if (local_freq.count(i)) {
                f = local_freq[i];
            }
            fp_t fsum = 0.0;
            MPI_Allreduce(&f, &fsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sample_hw_freq[i] = fsum;
        }
        MPI_Allreduce(&local_errors, &sample_errors, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_statistical_shots, &sample_statistical_shots, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&hwsum, &sample_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sample_mean /= shots_this_round;
        if (first_iter) {
            sample_statistical_shots = shots_this_round;
        }
#ifdef STATBENCH_DEBUG
        if (sample_statistical_shots < prev_sample_ss && world_rank == 0) {
            std::cout << "ss < prev_ss: " << sample_statistical_shots << " < " << prev_sample_ss << "\n";
        }
#endif
        // Update hamming_weight_freq table
        for (auto pair : hamming_weight_freq) {
            // Apply prior probabilities to "fill up" other shots (the ss - shots_per_round shots).
            hamming_weight_freq[pair.first] *= sample_statistical_shots/prev_sample_ss;
        }
        fp_t true_sample_ss = 0.0;  // To account for floating point error.
        for (auto pair : sample_hw_freq) {
            if (pair.second == 0) {
                continue;
            }
            if (!hamming_weight_freq.count(pair.first)) {
                hamming_weight_freq[pair.first] = 0;
            }
            hamming_weight_freq[pair.first] += pair.second;
            true_sample_ss += hamming_weight_freq[pair.first];
        }
        // Record statistics and delete current decoder afterward.
        statres.n_logical_errors = 
            statres.n_logical_errors * sample_statistical_shots/prev_sample_ss + local_errors;
        statres.statistical_shots = true_sample_ss;
        delete curr_decoder;
        // Update noise and create a new decoder.
        prev_sample_ss = sample_statistical_shots;
#ifdef STATBENCH_DEBUG
        if (world_rank == 0) {
            std::cout << "\tphysical error rate: " << pc << "\n";
            std::cout << "\tmean hamming weight: " << sample_mean << "\n";
            std::cout << "\tstatistical shots: " << sample_statistical_shots << "\n";
            std::cout << "\tdistribution:\n";
            fp_t psum = 0.0;
            fp_t shotsum = 0.0;
            for (auto pair : hamming_weight_freq) {
                if (pair.second) {
                    std::cout << "\t\t" << pair.first << "\t" 
                        << (pair.second / statres.statistical_shots) 
                        << "(" << pair.second << " of " << statres.statistical_shots << ")\n";
                    psum += pair.second / statres.statistical_shots;
                    shotsum += pair.second;
                }
            }
            std::cout << "\tprob sum: " << psum << "\n";
            std::cout << "\tshot sum: " << shotsum << " (should be " << statres.statistical_shots << ")\n";
        }
#endif
        pc = pow(M_E, log(pc) + log_p_update);
        curr_decoder = mkdec(pc);
        shots -= shots_this_round;
        first_iter = false;
    }
    delete curr_decoder;
    for (auto pair : hamming_weight_freq) {
        statres.hamming_weight_dist[pair.first] = pair.second / statres.statistical_shots;
    }

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
    fp_t reset_leak_mean,
    fp_t round_leak_stddev,
    fp_t clifford_leak_stddev,
    fp_t reset_leak_stddev)
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
    params.after_reset_leakage_probability = __CHS(error_mean, pauliplus_error_mean, reset_leak_mean);
    
    params.before_round_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, round_leak_stddev);
    params.after_clifford_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, clifford_leak_stddev);
    params.after_reset_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, reset_leak_stddev);

    params.both_stabilizers = both_stabilizers;
    params.swap_lru = other_flags & BC_FLAG_SWAP_LRU_V1;
    params.swap_lru_with_no_swap = other_flags & BC_FLAG_SWAP_LRU_V2;
    params.initial_state_is_basis_1 = other_flags & BC_FLAG_INVERT_STATE;

    stim::Circuit circ = generate_surface_code_circuit(params).circuit;
    return circ;
}

}  // qrc
