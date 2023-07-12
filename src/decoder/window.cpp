/*
 *  author: Suhas Vittal
 *  date:   7 July 2023
 * */

#include "decoder/window.h"

namespace qontra {
namespace decoder {

Decoder::result_t
WindowDecoder::decode_error(const syndrome_t& syndrome) {
    uint rounds = (circuit.count_detectors() / detectors_per_round) - 1;

    const uint det = base_decoder->get_circuit().count_detectors();
    const uint obs = base_decoder->get_circuit().count_observables();

    syndrome_t running_syndrome(syndrome);

    fp_t exec_time = 0;
    std::vector<Decoder::assign_t> error_assignments;
    stim::simd_bits corr(obs);
    for (uint r = 0; r < rounds; r += commit_window) {
        bool is_last_call = r + commit_window >= rounds;

        uint offset = r*detectors_per_round;
        uint upp = det - (is_last_call ? 0 : detectors_per_round);

        syndrome_t rxs(det);
        rxs.clear();
        for (uint i = 0; i < upp; i++) {
            rxs[i] = running_syndrome[i+offset];
        }
        Decoder::result_t rxres;
        if (is_last_call)   rxres = decode_final_error(rxs);
        else if (r == 0)    rxres = decode_first_error(rxs);
        else                rxres = decode_middle_error(rxs);

        for (auto aa : rxres.error_assignments) {
            uint x = std::get<0>(aa), y = std::get<1>(aa);
            running_syndrome[x + offset] ^= 1;
            running_syndrome[y + offset] ^= 1;

            error_assignments.push_back(
                    std::make_tuple(x+offset, y+offset, std::get<2>(aa)));
        }
        corr ^= rxres.corr;
        exec_time += rxres.exec_time;
    }
    
    Decoder::result_t final_res = {
        exec_time,
        corr,
        is_error(corr, syndrome),
        error_assignments
    };
    return final_res;
}

Decoder::result_t
WindowDecoder::decode_first_error(const syndrome_t& syndrome) {
    syndrome_t base_syndrome(base_decoder->get_circuit().count_detectors());
    base_syndrome.clear();
    base_syndrome |= syndrome;
    return retrieve_result_from_base(
                base_syndrome, 0, commit_window*detectors_per_round);
}

Decoder::result_t
WindowDecoder::decode_middle_error(const syndrome_t& syndrome) {
    syndrome_t base_syndrome(base_decoder->get_circuit().count_detectors());
    base_syndrome.clear();
    for (uint k : get_nonzero_detectors(syndrome)) {
        base_syndrome[detectors_per_round + k] = syndrome[k];
    }
    return retrieve_result_from_base(
                            base_syndrome, 
                            detectors_per_round,
                            (commit_window+1) * detectors_per_round);
}

Decoder::result_t
WindowDecoder::decode_final_error(const syndrome_t& syndrome) {
    syndrome_t base_syndrome(base_decoder->get_circuit().count_detectors());
    base_syndrome.clear();
    for (uint k : get_nonzero_detectors(syndrome)) {
        base_syndrome[detectors_per_round + k] = syndrome[k];
    }
    return retrieve_result_from_base(
                            base_syndrome, 
                            detectors_per_round,
                            std::numeric_limits<uint>::max());
}

Decoder::result_t
WindowDecoder::retrieve_result_from_base(
                            const syndrome_t& base_syndrome,
                            uint min_detector,
                            uint max_detector) 
{
    auto res = base_decoder->decode_error(base_syndrome);

    std::vector<Decoder::assign_t> error_assignments;
    stim::simd_bits corr(base_decoder->get_circuit().count_observables());
    corr.clear();

    // We only care about assignments in the first round.
    for (auto aa : res.error_assignments) {
        uint di = std::get<0>(aa);
        if (di < min_detector || di >= max_detector) {
            continue;
        }
        uint dj = std::get<1>(aa);
        auto local_corr = std::get<2>(aa);
        error_assignments.push_back(aa);
        corr ^= local_corr;
    }
    Decoder::result_t final_res = {
        res.exec_time,
        corr,
        false,  // No one cares about this value for sliding window decoding.
        error_assignments
    };
    return final_res;
}

}   // decoder
}   // qontra
