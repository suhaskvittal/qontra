/*
 *  author: Suhas Vittal
 *  date:   16 July 2023
 * */

#ifndef MATCHING_SYNDROME_PREDICTOR_h
#define MATCHING_SYNDROME_PREDICTOR_h

#include "decoder/decoder.h"
#include "sim/components/syndrome_predictor.h"

namespace qontra {
namespace sim {

// The following predictor leverages the insight that time-like
// errors will occur if there is a "long" assignment of one syndrome
// bit to another or the boundary. Thus, this predictor calls a 1-cycle decoder
// such as Astrea to quickly decode a small region of syndrome bits
// (does not care about other syndrome bits -- no need for perfect
// accuracy here). Then, we check if any assignments are unsatisfactory,
// or overly long over some threshold. If an assignment travels through
// the boundary, it is also considered unsatisfactory, but it may be
// doubly unsatisfactory (unsatisfactory for both syndrome bits) or
// singly unsatisfactory.
//
// If the Hamming weight of a region is too high, the decoder will default
// to no prediction.
//

class MatchingSyndromePredictor : public SyndromePredictor {
public:
    MatchingSyndromePredictor(
            uint round_depth,
            uint detectors_per_round,
            decoder::Decoder* dec,
            uint threshold,
            uint hamming_weight_limit=6)   // 10 is the default for Astrea
        :SyndromePredictor(round_depth, detectors_per_round),
        backing_decoder(dec),
        threshold(threshold),
        hamming_weight_limit(hamming_weight_limit)
    {}

    SyndromePredictor::pred_result_t    predict(stim::simd_bit_table&) override;
private:
    decoder::Decoder*   backing_decoder;
    uint                threshold;
    uint                hamming_weight_limit;
};

}   // sim
}   // qontra

#endif  // MATCHING_SYNDROME_PREDICTOR_h
