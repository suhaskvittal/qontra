/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODER_MWPM_h
#define DECODER_MWPM_h

#include "decoder/decoder.h"
#include "graph/decoding_graph.h"
#include "graph/graph.h"

#include <PerfectMatching.h>

#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <time.h>

namespace qontra {

typedef uint32_t    wgt_t;

const fp_t MWPM_INTEGER_SCALE = 1000.0;
#define MWPM_TO_INT(x)  ((wgt_t) (MWPM_INTEGER_SCALE * (x)))

class MWPMDecoder : public Decoder {
public:
    MWPMDecoder(DetailedStimCircuit circ)
        :Decoder(circ)
    {}

    std::string name(void) override { return "MWPMDecoder"; }

    Decoder::result_t       decode_error(const syndrome_t&) override;
private:
    //
    // Extra functionality for restriction decoder:
    //
    friend class RestrictionDecoder;
    std::map<std::pair<graph::decoding::vertex_t*, graph::decoding::vertex_t*>, fp_t> 
        override_weights;
};

}   // qontra

#endif  // DECODER_MWPM_h
