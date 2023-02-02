/*
 *  author: Suhas Vittal
 *  date:   9 August 2022
 * */

#include "benchmark.h"

namespace qrc {

#define STATBENCH_DEBUG
//#define STATBENCH_DEBUG2

/* Helper Functions */
inline fp_t __lognCk(fp_t n, fp_t k) {
    return lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1);
}

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
        stim::simd_bit_table sample_buffer = 
            stim::detector_samples(decoder_p->circuit, shots_this_round,
                    false, true, rng);
        // Last part of samples is the actual observable.
        // We are trying to match that.
        sample_buffer = sample_buffer.transposed();
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
                n_logical_errors += res.is_logical_error ? 1 : 0;
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

StatisticalResult
b_statistical_ler(dgf_t& mkdec, uint code_dist, fp_t p, uint64_t shots, std::mt19937_64& rng, uint64_t update_rate) {
    const fp_t max_uncorrectable_nlogprob = 11.512925465;  // Represents uncorrectable probability 1e-3;
    
    Decoder * curr_decoder = mkdec(p);
    // Compute mean physical error rate from Decoding Graph.
    fp_t mean_flip_prob = 0.0;
    uint64_t n_error_sources = 0;
    for (auto v : curr_decoder->graph.vertices()) {
        for (auto w : curr_decoder->graph.adjacency_list(v)) {
            auto e = curr_decoder->graph.get_edge(v, w);
            mean_flip_prob += e->error_probability;
            n_error_sources += 1;
        }
    }
    mean_flip_prob /= (fp_t)n_error_sources;
    const fp_t r = p / mean_flip_prob;
    const fp_t min_uncorrectable_nlogprob = error_prob_to_uncorrectable_nlogprob(code_dist, mean_flip_prob);

    // How many updates (multiplicative) we need to reach the max uncorrectable probability.
    const fp_t pucupdate = 
        ((fp_t)update_rate) * (min_uncorrectable_nlogprob - max_uncorrectable_nlogprob) / ((fp_t)shots);
    const fp_t ehw = n_error_sources * mean_flip_prob;
    const fp_t logehwprob = __lognCk(n_error_sources, ehw)
                            + ehw * log(mean_flip_prob)
                            + (n_error_sources - ehw) * log(1.0-mean_flip_prob);

#ifdef STATBENCH_DEBUG
    std::cout << "mean flip prob: " << mean_flip_prob << "\n";
    std::cout << "min uncorrectable nlogprob = " << min_uncorrectable_nlogprob << "\n";
    std::cout << "puc update = " << pucupdate << "\n";
    std::cout << "ratio: " << r << "\n";
    std::cout << "error sources: " << n_error_sources << "\n";
    std::cout << "expected hamming weight: " << ehw << " logprob = " << logehwprob << "\n";
#endif

    StatisticalResult statres;
    statres.true_shots = shots;

    std::map<uint, fp_t> hamming_weight_freq;

    fp_t pc = p;    // Current physical error rate
    fp_t nlogpuc = min_uncorrectable_nlogprob;  // Current uncorrectable probability
    fp_t prev_ss = shots < update_rate ? shots : update_rate;

    bool first_iter = true;  // We need to adjust any inaccuracies in our model in the first iteration only.
    while (shots > 0) {
#ifdef STATBENCH_DEBUG
        std::cout << "shots left: " << shots << "\n";
#endif
        const uint64_t shots_this_round = shots < update_rate ? shots : update_rate;
        // Benchmark decoder
        uint64_t s = shots_this_round;

        std::map<uint, fp_t> local_freq;
        uint64_t local_errors = 0;

        fp_t ss = 0.0;
        while (s > 0) {
            uint64_t shots_this_batch = s < 100'000 ? s : 100'000;
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
                if (!local_freq.count(hw)) {
                    local_freq[hw] = 0;
                }
                local_freq[hw]++;
                if (hw) {
                    auto res = curr_decoder->decode_error(syndrome);
                    local_errors += res.is_logical_error;
                    statres.mean_execution_time += res.execution_time;
                    if (res.execution_time > statres.max_execution_time) {
                        statres.max_execution_time = res.execution_time;
                    }
                }
                // Compute statistical shots for this batch.
                // Probability of achieving a certain Hamming weight follows Bin(n_error_sources, mean_flip_prob).
                const fp_t hwdiv2 = hw * 0.5;
                fp_t loghwprob = __lognCk(n_error_sources, hwdiv2)
                                + hwdiv2 * log(mean_flip_prob)
                                + ( ((fp_t)n_error_sources) - hwdiv2 ) * log(1.0-mean_flip_prob);
                fp_t logss = logehwprob - loghwprob;
                ss += pow(M_E, logss);
            }
            
            s -= shots_this_batch;
        }
        if (first_iter) {
            ss = shots_this_round;
        }
        // Update hamming_weight_freq table
        for (auto pair : hamming_weight_freq) {
            // Apply prior probabilities to "fill up" other shots (the ss - shots_per_round shots).
            hamming_weight_freq[pair.first] = pair.second * (ss - shots_this_round)/prev_ss;
        }
        for (auto pair : local_freq) {
            if (!hamming_weight_freq.count(pair.first)) {
                hamming_weight_freq[pair.first] = 0;
            }
            hamming_weight_freq[pair.first] += pair.second;
        }
        // Record statistics and delete current decoder afterward.
        statres.n_logical_errors = statres.n_logical_errors * (ss - shots_this_round)/prev_ss + local_errors;
        statres.statistical_shots += ss;
        delete curr_decoder;
        // Update noise and create a new decoder.
        prev_ss = ss;
        nlogpuc -= pucupdate;
        pc = uncorrectable_nlogprob_to_error_prob(code_dist, nlogpuc) * r;
#ifdef STATBENCH_DEBUG
        std::cout << "\tphysical error rate: " << pc << "\n";
        std::cout << "\tstatistical shots: " << ss << "\n";
#endif
        curr_decoder = mkdec(pc);
        shots -= shots_this_round;
        first_iter = false;
    }
    statres.mean_execution_time /= statres.statistical_shots;
#ifdef STATBENCH_DEBUG
    std::cout << "Distribution:\n";
#endif
    for (auto pair : hamming_weight_freq) {
        statres.hamming_weight_dist[pair.first] = pair.second / statres.statistical_shots;
#ifdef STATBENCH_DEBUG
        if (pair.second) {
            std::cout << "prob(" << pair.first << ") = " << statres.hamming_weight_dist[pair.first] << "\n";
        }
#endif
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
    uint8_t swap_lru,
    uint rounds,
    fp_t clevel_error_mean,
    fp_t clevel_error_stddev,
    fp_t pauliplus_error_mean,
    fp_t pauliplus_error_stddev,
    fp_t round_dp_mean,
    fp_t clifford_dp_mean,
    fp_t reset_flip_mean,
    fp_t meas_flip_mean,
    fp_t round_dp_stddev,
    fp_t clifford_dp_stddev,
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
    params.after_clifford_depolarization = __CHS(error_mean, clevel_error_mean, clifford_dp_mean);
    params.after_reset_flip_probability = __CHS(error_mean, clevel_error_mean, reset_flip_mean);
    params.before_measure_flip_probability = __CHS(error_mean, clevel_error_mean, meas_flip_mean);

    params.before_round_data_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, round_dp_stddev);
    params.after_clifford_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, clifford_dp_stddev);
    params.after_reset_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, reset_flip_stddev);
    params.before_measure_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, meas_flip_stddev);

    params.before_round_leakage_probability = __CHS(error_mean, pauliplus_error_mean, round_leak_mean);
    params.after_clifford_leakage_probability = __CHS(error_mean, pauliplus_error_mean, clifford_leak_mean);
    params.after_reset_leakage_probability = __CHS(error_mean, pauliplus_error_mean, reset_leak_mean);
    
    params.before_round_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, round_leak_stddev);
    params.after_clifford_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, clifford_leak_stddev);
    params.after_reset_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, reset_leak_stddev);

    params.both_stabilizers = both_stabilizers;
    params.swap_lru = swap_lru & 0b01;
    params.swap_lru_with_no_swap = swap_lru & 0b10;

    stim::Circuit circ = generate_surface_code_circuit(params).circuit;
    return circ;
}

}  // qrc
