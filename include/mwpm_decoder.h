/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef MWPM_DECODER_h
#define MWPM_DECODER_h

#include "decoder.h"

#include <PerfectMatching.h>

#include <vector>
#include <map>
#include <chrono>
#include <set>
#include <utility>

#include <time.h>

#define MWPM_DECODER_NAME "MWPMDecoder"
#define MWPM_INTEGER_SCALE 10000.0

typedef qfp_t wgt_t;

namespace qrc {

class MWPMDecoder : public Decoder {
public:
    MWPMDecoder(const stim::Circuit&, uint max_detector=BOUNDARY_INDEX);

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;

    uint32_t longest_error_chain;
protected:
    std::vector<uint8_t> get_correction_from_matching(const std::map<uint, uint>&);

    PathTable path_table;
    uint max_detector;
};

}  // qrc

#endif
