/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 * */

#ifndef SYNDROME_PREDICTOR_h
#define SYNDROME_PREDICTOR_h

#include "defs.h"

#include <stim.h>

namespace qontra {
namespace sim {

// SyndromePredictor is a virtual class for predicting syndrome extraction
// measurements in the next round. The input size for the predictor is
// round_depth*detectors_per_round and the output is detectors_per_round, which
// returns the predicted detection events.

class SyndromePredictor {
public:
    SyndromePredictor(
            uint round_depth,
            uint detectors_per_round)
        :round_depth(round_depth),
        detectors_per_round(detectors_per_round),
        trial_mask(1),
        trial_mask_valid(false)
    {}

    typedef struct {
        stim::simd_bit_table    sig_pred;   // Whether to skip the next round.
        stim::simd_bit_table    val_pred;   // The predicted measurement in 
                                            // the next round.

        uint64_t    cycles;     // Additional cycles required to compute the
                                // prediction
    } pred_result_t;

    virtual pred_result_t   predict(stim::simd_bit_table&) =0;

    const uint round_depth;
    const uint detectors_per_round;

    stim::simd_bits     trial_mask; // Use this to prevent predictors from
                                    // unnecessarily predicting on trials that
                                    // will be tossed away. Useful for predictors
                                    // that evaluate on a trial-by-trial basis.
    bool                trial_mask_valid;
};

// The SimplePredictor only predicts all-zeros detection events.

class SimpleSyndromePredictor : public SyndromePredictor {
public:
    SimpleSyndromePredictor(uint round_depth, uint detectors_per_round)
        :SyndromePredictor(round_depth, detectors_per_round)
    {}

    SyndromePredictor::pred_result_t
    predict(stim::simd_bit_table& input) override {
        stim::simd_bit_table sig_pred(
                detectors_per_round, input.num_minor_bits_padded(), true);
        stim::simd_bit_table val_pred(
                detectors_per_round, input.num_minor_bits_padded());

        for (uint r = 0; r < round_depth; r++) {
            for (uint i = 0; i < detectors_per_round; i++) {
                sig_pred[i].for_each_word(
                        input[r*detectors_per_round + i],
                        [&] (auto& sig, auto& in)
                        {
                            sig &= ~in;
                        });
            }
        }
        return (SyndromePredictor::pred_result_t) {
            sig_pred,
            val_pred,
            0
        };
    }
};

}   // sim
}   // qontra

#endif  // SYNDROME_PREDICTOR_h
