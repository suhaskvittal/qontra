/*
 *  author: Suhas Vittal
 *  date:   22 March 2023
 * */

#ifndef WINDOW_WRAPPER_h
#define WINDOW_WRAPPER_h

#include "decoder.h"
#include "defs.h"
#include "mwpm_decoder.h"

#include <map>
#include <vector>

#include <stim.h>

namespace qrc {
namespace window {

class Wrapper : public MWPMDecoder {
public:
    Wrapper(const stim::Circuit&, Decoder*, uint total_rounds, uint window_size, uint detectors_per_round,
            uint max_correction_depth);

    std::string name(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;

    fp_t sum_hamming_weight_in_window;
    uint64_t total_downtime;
    uint64_t max_hamming_weight_in_window;
protected:
    void merge_matching(const std::map<uint, uint>&, uint offset, bool limit_to_first_round=true);

    Decoder * decoder;

    const uint total_rounds;
    const uint window_size;
    const uint detectors_per_round;
    const uint max_correction_depth;

    uint64_t shot_hw_sum;
    uint64_t n_merge_calls;

    std::vector<uint8_t> running_syndrome;
    std::map<uint, uint> running_matching;
};

}
}

#endif
