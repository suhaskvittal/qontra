/*
 *  author: Suhas Vittal
 *  date:   22 March 2023
 * */

#ifndef WINDOW_ADAPTIVE_h
#define WINDOW_ADAPTIVE_h

#include "decoder.h"
#include "defs.h"
#include "mwpm_decoder.h"
#include "window/wrapper.h"

#include <fstream>
#include <iostream>

namespace qrc {
namespace window {

class AdaptiveDecoder : public Wrapper {
public:
    AdaptiveDecoder(const stim::Circuit&, Decoder*, uint total_rounds, uint detectors_per_round, 
            uint max_correction_depth);

    std::string name(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;

    bool recording;
    std::ofstream record_out;

    fp_t sum_window_size;
    uint64_t largest_window_size;
protected:
    void decode_syndrome(const std::vector<uint8_t>&, uint offset, bool last_round=false);
    void record_window_data(uint window_size, uint hw);

    uint64_t shot_ws_sum;
    uint64_t window_size;

    std::vector<uint8_t> running_correction;
};

}
}

#endif
