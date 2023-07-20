/*
 *  author: Suhas Vittal
 *  date:   16 July 2023
 * */

#ifndef MULTI_LEVEL_SYNDROME_PREDICTOR_h
#define MULTI_LEVEL_SYNDROME_PREDICTOR_h

#include "sim/components/syndrome_predictor.h"

#include <array>

#include <stim.h>

namespace qontra {
namespace sim {

// This class is designed to use multiple predictors
// in a hierarchical fashion. That is,
//      predictor 0's decisions have priority
//      predictor 1's decisions have priority wherever predictor 0
//          did not signal any prediction.
//      etc.

template <int L=2>
class MultiLevelPredictor : public SyndromePredictor {
public:
    MultiLevelPredictor()
        :SyndromePredictor(0, 0),
        predictors()
    {
        predictors.fill(nullptr);
    }

    SyndromePredictor::pred_result_t
    predict(stim::simd_bit_table& input_data) override {
        if (predictors[0] == nullptr) {
            return (SyndromePredictor::pred_result_t) {
                stim::simd_bit_table(1, 1),
                stim::simd_bit_table(1, 1),
                0
            };
        }
        auto running_res = predictors[0]->predict(input_data);
        for (uint k = 1; k < L; k++) {
            if (predictors[k] == nullptr)   break;
            auto new_res = predictors[k]->predict(input_data);
            for (uint i = 0; 
                    i < running_res.sig_pred.num_major_bits_padded(); i++) 
            {
                running_res.sig_pred[i].for_each_word(
                        running_res.val_pred[i],
                        new_res.sig_pred[i],
                        new_res.val_pred[i],
                        [&] (auto& s1, auto& v1, auto& s2, auto& v2)
                        {
                            v1 = (v1 & s1) | (v2 & ~s1 & s2);
                            s1 |= s2;
                        });
            }
        }
        return running_res;
    }

    SyndromePredictor*& operator[](int index) {
        assert(index >= 0 && index < L);
        return predictors[index];
    }
private:
    std::array<SyndromePredictor*, L> predictors;
};

}   // sim
}   // qontra

#endif  // MULTI_LEVEL_SYNDROME_PREDICTOR_h
