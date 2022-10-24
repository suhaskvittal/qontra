/*
 *  author: Suhas Vittal
 *  date:   1 October 2022
 * */

#ifndef GULLIVER_TRIAL_DECODER_h
#define GULLIVER_TRIAL_DECODER_h

#include "defs.h"
#include "mwpm_decoder.h"

#include <deque>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

#include <time.h>

#define BDC_DEBUG

namespace qrc {
namespace gulliver {

class TrialDecoder : public MWPMDecoder {
public:
    TrialDecoder(const stim::Circuit&);
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
private:
    std::map<uint, uint> hungarian(const std::vector<uint>&);
    std::map<uint, uint> blossom(const std::vector<uint>&);
};

}   // gulliver
}   // qrc

#endif
