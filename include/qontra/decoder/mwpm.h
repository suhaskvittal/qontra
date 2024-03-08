/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODER_MWPM_h
#define DECODER_MWPM_h

#include "qontra/decoder/matching_base.h"

namespace qontra {

class MWPMDecoder : public MatchingBase {
public:
    MWPMDecoder(const DetailedStimCircuit&);

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
};

}   // qontra

#endif  // DECODER_MWPM_h
