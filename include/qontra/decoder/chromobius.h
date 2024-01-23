/*
 *  author: Suhas Vittal
 *  date:   23 December 2023
 *
 *  Wrapper for Gidney and Jones's Chromobius decoder.
 * */

#ifndef QONTRA_CHROMOBIUS_h
#define QONTRA_CHROMOBIUS_h

#include "qontra/decoder.h"

#include <chromobius/decode/decoder.h>

namespace qontra {

chromobius::Decoder init_chromobius(stim::Circuit);

class Chromobius : public Decoder {
public:
    Chromobius(DetailedStimCircuit circuit)
        :Decoder(circuit, graph::DecodingGraph::Mode::DO_NOT_BUILD),
        backing_decoder(init_chromobius(circuit))
    {}
        
    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
private:
    chromobius::Decoder backing_decoder;
};

}   // qontra

#include "chromobius.inl"

#endif  // QONTRA_CHROMOBIUS_h
