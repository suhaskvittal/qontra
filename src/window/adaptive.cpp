/*
 *  author: Suhas Vittal
 *  date:   22 March 2023
 * */

#include "window/adaptive.h"

namespace qrc {
namespace window {

AdaptiveDecoder::AdaptiveDecoder(const stim::Circuit& circuit, Decoder * backing_decoder, uint total_rounds,
                                uint detectors_per_round, uint max_correction_depth)
    :Wrapper(circuit, backing_decoder, total_rounds, total_rounds, detectors_per_round, max_correction_depth),
    recording(false),
    sum_window_size(0),
    largest_window_size(0),
    shot_ws_sum(0),
    window_size(0),
    running_correction(circuit.count_observables())
{}

std::string
AdaptiveDecoder::name() {
    return "Adaptive(" + decoder->name() + ")";
}

DecoderShotResult
AdaptiveDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    // Decode on different partitions of the syndrome.
    running_correction = std::vector<uint8_t>(circuit.count_observables(), 0);
    running_syndrome = std::vector<uint8_t>(syndrome);
    running_matching.clear();

    std::vector<uint8_t> acc_syndrome;

    uint lower_d = 0;
    uint upper_d = detectors_per_round;
    uint lim = total_rounds * detectors_per_round;
    // State variables
    bool sig_seeking = false;
    uint rounds_since_hit = 0;

    window_size = 0;

    shot_ws_sum = 0;
    shot_hw_sum = 0;
    n_merge_calls = 0;

    while (upper_d < lim) {
        uint64_t hw = 0;
        uint k = upper_d - detectors_per_round;
        for (uint i = k; i < upper_d; i++) {
            acc_syndrome.push_back(running_syndrome[i]);
            hw += running_syndrome[i];
        }

        if (hw == 0) {
            if (sig_seeking) {
                rounds_since_hit++;
            } else {
                acc_syndrome.clear();
                lower_d += detectors_per_round;
                upper_d += detectors_per_round;
                window_size = 0;
                continue;
            }
        } else {
            sig_seeking = true;
            rounds_since_hit = 0;
        }
        window_size++;
        
        if (rounds_since_hit == max_correction_depth) {
            // Fit the decoder into the decoder parameters.
            decode_syndrome(acc_syndrome, lower_d);
            acc_syndrome.clear();
            sig_seeking = false;

            rounds_since_hit = 0;
            window_size = 0;
            lower_d = upper_d;
        }
        upper_d += detectors_per_round;
    }
    // Decode final round.
    // Push the remaining bits onto acc_syndrome.
    window_size++;
    uint k = upper_d - detectors_per_round;
    for (uint i = k; i < circuit.count_detectors(); i++) {
        acc_syndrome.push_back(running_syndrome[i]);
    }
    uint hw = std::accumulate(acc_syndrome.begin(), acc_syndrome.end(), 0);
    if (hw) {
        decode_syndrome(acc_syndrome, 0, true);
    }
    auto correction = get_correction_from_matching(running_matching);
    bool is_error = is_logical_error(correction, syndrome, 
                        circuit.count_detectors(), circuit.count_observables());
    if (is_error) {
        auto mwpm_res = MWPMDecoder::decode_error(syndrome);
        if (!mwpm_res.is_logical_error) {
            std::cout << "=================================\n";
            // Compare matchings
            std::cout << "window matching:\n";
            for (auto pair : running_matching) {
                int r1 = pair.first == BOUNDARY_INDEX ? -1 : pair.first/detectors_per_round;
                int r2 = pair.second == BOUNDARY_INDEX ? -1 : pair.second/detectors_per_round;
                std::cout << "\t" << pair.first << "(" << r1 << ")"
                            << " --> " << pair.second << "(" << r2 << ")" << "\n";
            }

            std::cout << "mwpm matching:\n";
            for (auto pair : mwpm_res.matching) {
                int r1 = pair.first == BOUNDARY_INDEX ? -1 : pair.first/detectors_per_round;
                int r2 = pair.second == BOUNDARY_INDEX ? -1 : pair.second/detectors_per_round;
                std::cout << "\t" << pair.first << "(" << r1 << ")"
                            << " --> " << pair.second << "(" << r2 << ")" << "\n";
            }
        }
    }

    if (n_merge_calls) {
        sum_hamming_weight_in_window += ((fp_t)shot_hw_sum) / ((fp_t) n_merge_calls);
        sum_window_size += ((fp_t)shot_ws_sum) / ((fp_t) n_merge_calls);
    }

    auto final_res = (DecoderShotResult) {
        0.0,
        0.0,
        is_error,
        correction,
        running_matching,
    };
    return final_res;
}

void
AdaptiveDecoder::decode_syndrome(const std::vector<uint8_t>& acc_syndrome, uint offset, bool last_round) {
    uint64_t hw = std::accumulate(acc_syndrome.begin(), acc_syndrome.end(), 0);
    shot_hw_sum += hw;
    shot_ws_sum += window_size;

    if (window_size > largest_window_size) {
        largest_window_size = window_size;
    }

    if (hw == 0) {
        return;
    } else if (hw > max_hamming_weight_in_window) {
        max_hamming_weight_in_window = hw;
    }
    record_window_data(window_size, hw);

    std::vector<uint8_t> dec_syndrome(detectors_per_round, 0);
    uint syndrome_size = decoder->circuit.count_detectors() + decoder->circuit.count_observables();
    if (last_round) {
        while (dec_syndrome.size() < decoder->circuit.count_detectors() - acc_syndrome.size()) {
            dec_syndrome.push_back(0);
        }
        for (uint8_t b : acc_syndrome) {
            dec_syndrome.push_back(b);
        }
        for (uint i = 0; i < decoder->circuit.count_observables(); i++) {
            dec_syndrome.push_back(0);
        }
    } else {
        for (uint8_t b : acc_syndrome) {
            dec_syndrome.push_back(b);
        }
        while (dec_syndrome.size() < syndrome_size) {
            dec_syndrome.push_back(0);
        }
    }
    auto res = decoder->decode_error(dec_syndrome);
    merge_matching(res.matching, offset, false);
}

void
AdaptiveDecoder::record_window_data(uint ws, uint hw) {
    if (!recording) {
        return;
    }
    out << ws << "," << hw << "\n";
}

}   // window
}   // qrc
