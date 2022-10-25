/*
 *  author: Suhas Vittal
 *  date:   9 August 2022
 * */

#include "benchmark.h"

namespace qrc {

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
b_decoder_ler(Decoder * decoder_p, uint32_t shots, std::mt19937_64& rng,
        bool save_per_shot_data) 
{
    // Clear stats.
    decoder_p->clear_stats();
    // Declare statistics
    uint32_t array_size = save_per_shot_data ? shots : 1;

    uint32_t total_shots = shots;
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

std::vector<uint8_t> _to_vector(const stim::simd_bits_range_ref& array,
        uint n_detectors, uint n_observables) 
{
    std::vector<uint8_t> syndrome(n_detectors+n_observables);
    for (uint i = 0; i < n_detectors + n_observables; i++) {
        syndrome[i] = array[i];
    }
    return syndrome;
}

}  // qrc
