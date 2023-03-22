/*
 *  author: Suhas Vittal
 *  date:   22 March 2023
 * */

#include "window/wrapper.h"

namespace qrc {
namespace window {

Wrapper::Wrapper(const stim::Circuit& circuit, Decoder * backing_decoder, uint total_rounds, uint window_size,
                    uint detectors_per_round)
    :MWPMDecoder(circuit),
    decoder(backing_decoder),
    total_rounds(total_rounds),
    window_size(window_size),
    detectors_per_round(detectors_per_round),
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

    while (upper_d < lim) {
        std::vector<uint8_t> sub_syndrome((window_size+1)*detectors_per_round+circuit.count_observables(), 0);
        for (uint i = lower_d; i < upper_d; i++) {
            sub_syndrome[i-lower_d] = running_syndrome[i];
        }
        auto res = decoder->decode_error(sub_syndrome);
        merge_matching(res.matching, lower_d);

        lower_d += detectors_per_round;
        upper_d += detectors_per_round;
    }
    // Last window is special as this is the last set.
    std::vector<uint8_t> sub_syndrome(running_syndrome.begin()+lower_d,
                                        running_syndrome.end()-circuit.count_observables());
    auto res = decoder->decode_error(sub_syndrome);
    merge_matching(res.matching, lower_d, false);

    auto correction = get_correction_from_matching(running_matching);
    bool is_error = is_logical_error(correction, syndrome, circuit.count_detectors(), circuit.count_observables());

    if (is_error) {
        std::cout << "=================================\n";
        auto mwpm_res = MWPMDecoder::decode_error(syndrome);
        // Compare matchings
        std::cout << "window matching:\n";
        for (auto pair : running_matching) {
            std::cout << "\t" << pair.first << " --> " << pair.second << "\n";
        }

        std::cout << "mwpm matching:\n";
        for (auto pair : mwpm_res.matching) {
            std::cout << "\t" << pair.first << " --> " << pair.second << "\n";
        }
    }

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
    for (auto pair : src) {
        // Only match detectors in the first round of the window
        if (pair.first > pair.second) {
            continue;
        }
        if (limit_to_first_round && pair.first > detectors_per_round) {
            continue;
        }
        uint adj_d1 = pair.first == BOUNDARY_INDEX ? BOUNDARY_INDEX : offset + pair.first;
        uint adj_d2 = pair.second == BOUNDARY_INDEX ? BOUNDARY_INDEX : offset + pair.second;
        if (adj_d1 == BOUNDARY_INDEX && running_matching.count(BOUNDARY_INDEX)) {
            // Merge the matching, so to speak.
            adj_d1 = running_matching[adj_d1];
            running_matching.erase(BOUNDARY_INDEX);
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
