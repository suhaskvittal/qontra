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
    Wrapper(const stim::Circuit&, Decoder*, uint total_rounds, uint window_size, uint detectors_per_round);

    std::string name(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
protected:
    void merge_matching(const std::map<uint, uint>&, uint offset, bool limit_to_first_round=true);

    Decoder * decoder;

    const uint total_rounds;
    const uint window_size;
    const uint detectors_per_round;

    std::vector<uint8_t> running_syndrome;
    std::map<uint, uint> running_matching;
};

}
}

#endif
