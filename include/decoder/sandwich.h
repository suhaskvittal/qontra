/*
 *  author: Suhas Vittal
 *  date:   19 September 2023
 * */

#ifndef DECODER_SANDWICH_h
#define DECODER_SANDWICH_h

#include "decoder/decoder.h"

namespace qontra {

class SandwichDecoder : public Decoder {
public:
    SandwichDecoder(const stim::Circuit& target_circuit,
                        Decoder* base,
                        uint commit_window,
                        uint buffer_window,
                        uint detectors_per_round)
        :Decoder(target_circuit, graph::DecodingGraph::Mode::LOW_MEMORY),
        base_decoder(base),
        commit_window(commit_window),
        buffer_window(buffer_window),
        detectors_per_round(detectors_per_round)
    {}

protected:

};

}   // qontra

#endif  // DECODER_SANDWICH_h
