/*
 *  author: Suhas Vittal
 *  date:   22 March 2023
 * */

#include "window/wrapper.h"

namespace qrc {
namespace window {

Wrapper::Wrapper(const stim::Circuit& circuit, Decoder * backing_decoder, uint total_rounds, uint window_size,
                    uint detectors_per_round, uint max_correction_depth)
    :MWPMDecoder(circuit),
    total_downtime(0),
    sum_hamming_weight_in_window(0),
    max_hamming_weight_in_window(0),
    decoder(backing_decoder),
    total_rounds(total_rounds),
    window_size(window_size),
    detectors_per_round(detectors_per_round),
    max_correction_depth(max_correction_depth),
    shot_hw_sum(0),
    n_merge_calls(0),
    running_syndrome(),
    running_matching()
{}

std::string
Wrapper::name() {
    return "Window(" + decoder->name() + ")";
}

DecoderShotResult
Wrapper::decode_error(const std::vector<uint8_t>& syndrome) {
    uint lower_d = 0;
    uint upper_d = detectors_per_round * window_size;
    uint lim = detectors_per_round * total_rounds;

    running_syndrome = std::vector<uint8_t>(syndrome);
    running_matching.clear();

    shot_hw_sum = 0;
    n_merge_calls = 0;

    while (upper_d < lim) {
        std::vector<uint8_t> sub_syndrome((window_size+2)*detectors_per_round+circuit.count_observables(), 0);
        uint64_t hw = 0;
        for (uint i = lower_d; i < upper_d; i++) {
            sub_syndrome[i-lower_d + detectors_per_round] = running_syndrome[i];
            hw += running_syndrome[i];
        }
        if (hw == 0) {
            total_downtime++;
            lower_d += detectors_per_round;
            upper_d += detectors_per_round;
            continue;
        }

        if (hw > max_hamming_weight_in_window) {
            max_hamming_weight_in_window = hw;
        }
        shot_hw_sum += hw;

        auto res = decoder->decode_error(sub_syndrome);
        merge_matching(res.matching, lower_d);
        lower_d += detectors_per_round;
        upper_d += detectors_per_round;
    }
    // Last window is special as this is the last set.
    std::vector<uint8_t> sub_syndrome(detectors_per_round, 0);
    uint64_t hw = 0;
    for (uint i = lower_d; i < running_syndrome.size() - circuit.count_observables(); i++) {
        sub_syndrome.push_back(running_syndrome[i]);
        hw += running_syndrome[i];
    }
    if (hw) {
        auto res = decoder->decode_error(sub_syndrome);
        merge_matching(res.matching, lower_d, false);
        if (hw > max_hamming_weight_in_window) {
            max_hamming_weight_in_window = hw;
        }
        shot_hw_sum += hw;
    } else {
        total_downtime++;
    }

    if (n_merge_calls) {
        sum_hamming_weight_in_window += ((fp_t)shot_hw_sum) / ((fp_t)n_merge_calls);
    }

    auto correction = get_correction_from_matching(running_matching);
    bool is_error = is_logical_error(correction, syndrome, circuit.count_detectors(), circuit.count_observables());

    DecoderShotResult final_result = {
        0.0, // Todo
        0.0, // Todo
        is_error,
        correction,
        running_matching
    };
    return final_result;
}

void
Wrapper::merge_matching(const std::map<uint, uint>& src, uint offset, bool limit_to_first_round) {
    n_merge_calls++;

    for (auto pair : src) {
        // Only match detectors in the first round of the window
        if (pair.first > pair.second) {
            continue;
        }
        if (limit_to_first_round && pair.first > 2*detectors_per_round) {
            continue;
        }
        uint adj_d1 = offset + pair.first - detectors_per_round;
        uint adj_d2 = pair.second == BOUNDARY_INDEX ? BOUNDARY_INDEX : offset + pair.second - detectors_per_round;
        // Make sure these two are not too far apart.
        if (limit_to_first_round && pair.second != BOUNDARY_INDEX) {
            uint r1 = pair.first / detectors_per_round;
            uint r2 = pair.second / detectors_per_round;
            if (r2 - r1 > max_correction_depth) {
                // Set second detector to boundary.
                adj_d2 = BOUNDARY_INDEX;
            }
        }
        if (adj_d2 == BOUNDARY_INDEX && running_matching.count(BOUNDARY_INDEX)) {
            adj_d2 = running_matching[adj_d2];
            running_matching.erase(BOUNDARY_INDEX);
        }

        if (adj_d1 != BOUNDARY_INDEX) {
            running_syndrome[adj_d1] = 0;
        } 
        if (adj_d2 != BOUNDARY_INDEX) {
            running_syndrome[adj_d2] = 0;
        }
        running_matching[adj_d1] = adj_d2;
        running_matching[adj_d2] = adj_d1;
    }
}

}   // window
}   // qrc
