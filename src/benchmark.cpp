/*
 *  author: Suhas Vittal
 *  date:   9 August 2022
 * */

#include "benchmark.h"

namespace qrc {

/* Helper Functions */
fp_t __lognCk(uint64_t n, uint64_t k) {
    static std::map<std::pair<uint64_t, uint64_t>, uint64_t> table;
    if (k < n - k) {
        return __lognCk(n, n-k);
    }

    auto n_k = std::make_pair(n, k);
    if (table.count(n_k)) {
        return table[n_k];
    }

    fp_t fpn = (fp_t)n;
    fp_t fpk = (fp_t)k;
    // Use Stirling's approximation for factorials. 
    fp_t top = log(sqrt(2*M_PI*fpn)) + fpn*(log(fpn) - 1);
    fp_t bot1 = log(sqrt(2*M_PI*(fpn-fpk))) + (fpn-fpk)*(log(fpn-fpk) - 1);
    fp_t bot2 = log(sqrt(2*M_PI*fpk)) + fpk*(log(fpk) - 1);
    fp_t s = top - bot1 - bot2;
    table[n_k] = s;
    return s;
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
}

uint64_t
b_statistical_ler(dgf_t& mkdec, uint code_dist, fp_t p, uint64_t shots, std::mt19937_64& rng, uint64_t update_rate) {
    const fp_t max_uncorrectable_nlogprob = 3;  // Represents uncorrectable probability 1e-3;
    
    Decoder * curr_decoder = mkdec(p);
    // Compute mean physical error rate from Decoding Graph.
    fp_t mean_flip_prob = 0.0;
    uint64_t n_error_sources = 0;
    for (auto v : curr_decoder->graph.vertices()) {
        for (auto w : curr_decoder->graph.adjacency_list(v)) {
            auto e = curr_decoder->graph.get_edge(v, w);
            mean_flip_prob += e.error_probability;
            n_error_sources += 1;
        }
    }
    mean_flip_prob /= (fp_t)n_error_sources;
    const fp_t r = p / mean_flip_prob;
    const fp_t min_uncorrectable_nlogprob = error_prob_to_uncorrectable_nlogprob(mean_flip_prob);

    // How many updates (multiplicative) we need to reach the max uncorrectable probability.
    const fp_t pucupdate = 
        ((fp_t)update_rate) * (min_correctable_nlogprob - max_uncorrectable_nlogprob) / ((fp_t)shots);

    StatisticalResults statres;
    statres.true_shots = shots;

    fp_t pc = p;    // Current physical error rate
    fp_t nlogpuc = min_uncorrectable_nlogprob;  // Current uncorrectable probability
    while (shots > 0) {
        uint64_t shots_this_round = shots < update_rate ? shots : update_rate;
        // Benchmark decoder
        stim::simd_bit_table sample_buffer = 
            stim::detector_samples(curr_decoder->circuit, shots_this_round,
                    false, true, rng);
        sample_buffer = sample_buffer.transposed();
        fp_t mean_hamming_weight = 0;
        for (uint64_t i = 0; i < shots_this_round; i++) {
            std::vector<uint8_t> syndrome = _to_vector(sample_buffer[i], 
                                                    curr_decoder->circuit.count_detectors(),
                                                    curr_decoder->circuit.count_observables());
            uint hw = std::accumulate(syndrome.begin(),
                                    syndrome.begin() + curr_decoder->circuit.count_detectors(),
                                    0);
            if (hw & 0x1) {
                hw++;
            }
            mean_hamming_weight += hw;
            auto res = curr_decoder->decode_error(syndrome);
            statres.n_logical_errors += res.is_logical_error;
            statres.mean_execution_time += res.execution_time;
            if (res.execution_time > statres.max_execution_time) {
                statres.max_execution_time = res.execution_time;
            }
        }
        mean_hamming_weight /= (fp_t)shots_this_round;
        // Compute statistical shots for this batch.
        // Probability of achieving a certain Hamming weight follows Bin(n_error_sources, mean_flip_prob).
        fp_t loghwprob = __lognCk(n_error_sources, (uint64_t)mean_hamming_weight) 
                        + mean_hamming_weight * log(mean_flip_prob)
                        + ( ((fp_t)n_error_sources) - mean_hamming_weight ) * log(mean_flip_prob);
        fp_t logss = np.log( (fp_t) shots_this_round ) - loghwprob;
        fp_t ss = pow(M_E, logss); 
        // Record statistics and delete current decoder afterward.
        statres.statistical_shots += ss;
        statres.n_logical_errors += curr_decoder->n_logical_errors;
        statres.mean_execution_time += curr_decoder->mean_execution_time * ss;
        if (curr_decoder->max_execution_time > statres.max_execution_time) {
            statres.max_execution_time = curr_decoder->max_execution_time;
        }
        delete curr_decoder;
        // Update noise and create a new decoder.
        nlogpuc += pucupdate;
        pc = uncorrectable_nlogprob_to_error_prob(nlogpuc) * r;
        curr_decoder = mkdec(pc);
        shots -= shots_this_round;
    }
    statres.mean_execution_time /= statres.statistical_shots;

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

    std::cout << params.before_round_data_depolarization << "," << error_mean
                << "," << clevel_error_mean << "," << round_dp_mean << "\n";

    stim::Circuit circ = generate_surface_code_circuit(params).circuit;
    return circ;
}

}  // qrc
