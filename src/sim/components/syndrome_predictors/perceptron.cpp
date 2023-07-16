/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 * */

#include "sim/components/syndrome_predictors/perceptron.h"

namespace qontra {
namespace sim {

SyndromePredictor::pred_result_t
PerceptronPredictor::predict(stim::simd_bit_table& input) {
    stim::simd_bit_table sig_pred_tp(input.num_minor_bits_padded(),
                                    detectors_per_round);
    stim::simd_bit_table val_pred_tp(input.num_minor_bits_padded(),
                                    detectors_per_round);
    auto input_tp = input.transposed();
    
    const uint off = round_depth*detectors_per_round;
    for (uint64_t t = 0; t < input.num_minor_bits_padded(); t++) {
        if (trial_mask_valid && !trial_mask[t]) continue;
        for (uint i = 0; i < detectors_per_round; i++) {
            wgt_t pred = get_bias(i);

            for (uint r = 0; r < round_depth; r++) {
                for (uint j = 0; j < detectors_per_round; j++) {
                    int iv = input_tp[t][j] ? 1 : -1;
                    wgt_t& x = get(i, r, j);
                    pred += x*iv;
                }
            }
            bool last_round_was_zero = !input_tp[t][i];

            pred.collapse(threshold);
            if (pred.is_zero > 0) {
                if ((pred.is_same < 0 && last_round_was_zero)
                    || (pred.is_same > 0 && !last_round_was_zero))
                {
                    // This is a conflict.
                    sig_pred_tp[t][i] = 0;
                } else {
                    sig_pred_tp[t][i] = 1;
                    val_pred_tp[t][i] = 0;
                }
            } else if (pred.is_zero == 0) {
                sig_pred_tp[t][i] = pred.is_same != 0;
                val_pred_tp[t][i] = pred.is_same > 0 
                                        ? input[t][off+i] : !input[t][off+i];
            } else {
                if ((pred.is_same < 0 && !last_round_was_zero)
                    || (pred.is_same > 0 && last_round_was_zero))
                {
                    // This is a conflict
                    sig_pred_tp[t][i] = 0;
                } else {
                    sig_pred_tp[t][i] = 1;
                    val_pred_tp[t][i] = 1;
                }
            }
        }
    }

    return (SyndromePredictor::pred_result_t) {
        sig_pred_tp.transposed(),
        val_pred_tp.transposed(),
        0
    };
}

void
PerceptronPredictor::train(stim::simd_bit_table& data) {
    stim::simd_bit_table data_tp = data.transposed();

    const uint off1 = round_depth*detectors_per_round;
    const uint off2 = off1 + detectors_per_round;
    for (uint64_t t = 0; t < data.num_minor_bits_padded(); t++) {
        if (t % 10000 == 0)
        std::cout << "\t\tt = " << t << "\n";
        for (uint i = 0; i < detectors_per_round; i++) {
            // First, compute the labels.
            wgt_t y;
            y.is_zero = data_tp[t][off2+i] == 0 ? 1 : -1;
            y.is_same = data_tp[t][off1+i] == data_tp[t][off2+i] ? 1 : -1;

            wgt_t& bias = get_bias(i);
            bias += y;
            bias.clamp();
            // Do not update weights multiple times if they are used
            // multiple times.
            for (uint r = 0; r < round_depth; r++) {
                for (uint j = 0; j < detectors_per_round; j++) {
                    wgt_t& x = get(i, r, j);
                    int iv = data_tp[t][j] ? 1 : -1;
                    x += iv*y;
                    x.clamp();
                }
            }
        }
    }
}

PerceptronPredictor::wgt_t&
PerceptronPredictor::get_bias(uint i) {
    const uint off = detectors_per_round*round_depth + 1;
    if (mode == Mode::share) {
        return weights[0];
    } else if (mode == Mode::geom) {
        return weights[i*off];
    } else if (mode == Mode::local) {
        return weights[i*off];
    } else {
        return weights[i*off];
    }
}

PerceptronPredictor::wgt_t&
PerceptronPredictor::get(uint i, uint r, uint j) {
    auto acc = std::make_tuple(i, r, j);
    if (access_memoizer.count(acc)) {
        int index = access_memoizer[acc];
        if (index < 0) {
            index = -index;
            weights[index] *= ZERO_WEIGHT;
        }
        return weights[index];
    }

    const uint off = detectors_per_round*round_depth + 1;
    if (mode == Mode::share) {
        access_memoizer[acc] = 1 + r*detectors_per_round + j;
        return weights[1 + r*detectors_per_round + j];
    } else if (mode == Mode::geom) {
        auto v1 = lattice->get_qubit_at_meas_t(i);
        auto v2 = lattice->get_qubit_at_meas_t(j);
        auto data = lattice->get_path_data(v1, v2);
        access_memoizer[acc] = i*off + 1 + data.error_chain.size();
        return weights[i*off + 1 + data.error_chain.size()];
    } else if (mode == Mode::local) {
        auto v1 = lattice->get_qubit_at_meas_t(i);
        auto v2 = lattice->get_qubit_at_meas_t(j);
        auto data = lattice->get_path_data(v1, v2);
        if (data.error_chain.size() > locality_level) {
            weights[i*off + 1 + r*detectors_per_round + j] *= ZERO_WEIGHT;
            access_memoizer[acc] = -(i*off + 1 + r*detectors_per_round + j);
        } else {
            access_memoizer[acc] = i*off + 1 + r*detectors_per_round + j;
        }
        return weights[i*off + 1 + r*detectors_per_round + j];
    } else {
        access_memoizer[acc] = i*off + 1 + r*detectors_per_round + j;
        return weights[i*off + 1 + r*detectors_per_round + j];
    }
}

}   // sim
}   // qontra
