/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODER_MWPM_h
#define DECODER_MWPM_h

#include "qontra/decoder.h"
#include "qontra/graph/decoding_graph.h"

namespace qontra {

const fp_t MWPM_INTEGER_SCALE = 1000.0;
#define MWPM_TO_INT(x)  ((uint32_t) (MWPM_INTEGER_SCALE * (x)))

class MWPMDecoder : public Decoder {
public:
    MWPMDecoder(DetailedStimCircuit circ)
        :Decoder(circ)
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
};

}   // qontra

#endif  // DECODER_MWPM_h
