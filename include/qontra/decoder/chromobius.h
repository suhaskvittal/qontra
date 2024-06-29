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
    Chromobius(const DetailedStimCircuit&);
        
    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
private:
    chromobius::Decoder backing_decoder;
};

}   // qontra

#endif  // QONTRA_CHROMOBIUS_h
