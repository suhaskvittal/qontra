/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODER_MWPM_h
#define DECODER_MWPM_h

#include "decoder.h"
#include "decoding_graph.h"
#include "graph/dijkstra.h"

#include <PerfectMatching.h>

#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <time.h>

namespace qontra {
namespace decoder {

typedef uint32_t    wgt_t;

const fp_t MWPM_INTEGER_SCALE = 1000.0;
#define MWPM_TO_INT(x)  ((wgt_t) (MWPM_INTEGER_SCALE * (x)))

class MWPMDecoder : public Decoder {
public:
    MWPMDecoder(const stim::Circuit& circ)
        :Decoder(circ)
    {}

    std::string name(void) override { return "MWPMDecoder"; }
    bool        is_software(void) override { return true; }

    Decoder::result_t       decode_error(const vsyndrome_t&) override;
};

}   // decoder
}   // qontra

#endif  // DECODER_MWPM_h
