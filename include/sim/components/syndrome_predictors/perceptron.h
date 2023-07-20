/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 * */

#ifndef PERCEPTRON_SYNDROME_PREDICTOR_h
#define PERCEPTRON_SYNDROME_PREDICTOR_h

#include "graph/lattice_graph.h"
#include "sim/components/syndrome_predictor.h"

#include <set>
#include <vector>

#include <math.h>

namespace qontra {
namespace sim {

// PerceptronPredictor is a simple perceptron-based
// predictor, which predicts:
//  (1) All-zeros
//  (2) Same as the last round
//  (3) Neither/Conflict
// Thus, there are two perceptrons, which predict each
// value. If the perceptron's predictions conflict (i.e.
// all-zeros is predicted and same (not all-zeros) is also
// predicted), then the predictors defaults to "not taken".
//
// The perceptron predictor additionally has three modes:
//  (1) Share weights. There will be round_depth*detectors_per_round+1
//      weights used, and each prediction will use the same weights.
//  (2) Geometric weights. There are weights for each "ball" around a check.
//  (3) Local weights. There are only weights for nearby checks. locality_level
//      controls what is considered "nearby".
//  (4) Global weights. Each check has round_depth*detectors_per_round+1
//      weights. This is an idealized approach.

// These are upper bounds, not the exact number of weights needed.
#define __SQR(x)                ((x)*(x))
#define __MAX_WEIGHTS(r, dpr)   __SQR((r)*(dpr)+1)

class PerceptronPredictor : public SyndromePredictor {
public:
    enum class Mode { share, geom, local, global };

    PerceptronPredictor(uint round_depth,
                        uint detectors_per_round, 
                        uint threshold,
                        graph::LatticeGraph* lattice,
                        Mode mode,
                        uint locality_level=1)
        :SyndromePredictor(round_depth, detectors_per_round),
        weights(__MAX_WEIGHTS(round_depth, detectors_per_round), {0, 0}),
        threshold(threshold),
        lattice(lattice),
        mode(mode),
        locality_level(locality_level),
        access_memoizer(),
        predict_only_nontrivial_events(false)
    {}

    PerceptronPredictor(const PerceptronPredictor& other)
        :SyndromePredictor(other),
        weights(other.weights),
        threshold(other.threshold),
        lattice(other.lattice),
        mode(other.mode),
        locality_level(other.locality_level),
        access_memoizer(other.access_memoizer),
        predict_only_nontrivial_events(other.predict_only_nontrivial_events)
    {}

    SyndromePredictor::pred_result_t    predict(stim::simd_bit_table&) override;

    // The training data should have extra rows that contain events in
    // the next round.
    void    train(stim::simd_bit_table&);

    // This function
    void    train_using_error_model(const stim::Circuit& error_model, 
                                        uint64_t shots);

    void print_weights(void) {
        std::cout << "Weights:\n";
        for (auto& x : weights) {
            if (x.is_zero == 0 && x.is_same == 0)   continue;
            std::cout << "\t" << x.is_zero << ", " << x.is_same << "\n";
        }
    }

    bool    predict_only_nontrivial_events; // This restricts the Perceptron
                                            // to only training/predicting
                                            // on events where the most recent
                                            // bit is 1.
private:
    struct wgt_t {
        int     is_zero = 0;
        int     is_same = 0;

        void clamp(void) {
            std::vector<int*> data{&is_zero, &is_same};
            for (auto x : data) {
                if (abs(*x) > 128) *x = (*x < 0 ? -1 : 1) * 128;
            }
        }

        wgt_t operator+=(const wgt_t& other) {
            this->is_zero += other.is_zero;
            this->is_same += other.is_same;
            return *this;
        }

        wgt_t operator*=(const wgt_t& other) {
            this->is_zero *= other.is_zero;
            this->is_same *= other.is_same;
            return *this;
        }

        void collapse(uint threshold) {
            std::vector<int*> data{&is_zero, &is_same};
            for (auto x : data) {
                if (abs(*x) < threshold) *x = 0;
                else    *x = *x < 0 ? -1 : 1;
            }
        }
    };

    bool    is_trivial(stim::simd_bit_table&, uint64_t, uint);

    wgt_t&  get_bias(uint i);
    wgt_t&  get(uint i, uint r, uint j);

    std::vector<wgt_t>  weights;
    const uint threshold;

    graph::LatticeGraph* lattice;

    Mode mode;
    uint locality_level;

    typedef std::tuple<uint, uint, uint> access_t;
    std::map<access_t, int>    access_memoizer; // This will speed up weight
                                                // accesses if get is rather
                                                // slow.
        
    // Below are some operator overloads to make the code readable.
    friend wgt_t operator+(int lhs, wgt_t rhs) {
        return { rhs.is_zero + lhs, rhs.is_same + lhs};
    }

    friend wgt_t operator*(int lhs, wgt_t rhs) {
        return { rhs.is_zero * lhs, rhs.is_same * lhs};
    }

    friend wgt_t operator+(wgt_t lhs, int rhs) { return rhs + lhs; }
    friend wgt_t operator*(wgt_t lhs, int rhs) { return rhs * lhs; }

    friend wgt_t operator+(wgt_t lhs, wgt_t rhs) { return (lhs += rhs); }
    friend wgt_t operator*(wgt_t lhs, wgt_t rhs) { return (lhs *= rhs); }

    constexpr static wgt_t ZERO_WEIGHT = {0, 0};
};

}   // sim
}   // qontra

#endif  // PERCEPTRON_SYNDROME_PREDICTOR_h

