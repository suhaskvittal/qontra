/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#include "decoder/neural.h"

namespace qontra {

using namespace experiments;

static const uint64_t PER_HW_LIMIT = 100'000;

void
NeuralDecoder::train(uint64_t shots, bool verbose) {
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();
    arma::mat data_matrix(det, shots, arma::fill::value(-1.0));
    arma::mat y(obs, shots, arma::fill::zeros);

    uint64_t shots_elapsed = 0;
    callback_t cb;
    // We will train the model every batch.
    std::map<uint, uint64_t> hw_freq;
    cb.prologue = [&] (stim::simd_bits_range_ref x) {
        std::vector<uint> detectors = get_nonzero_detectors(x);
        const uint hw = detectors.size();
        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
            return;
        }
        if (hw_freq[hw]++ >= PER_HW_LIMIT)  return;
        for (uint d : detectors) {
            data_matrix(d, shots_elapsed) = 1;
        }
        for (uint i = 0; i < obs; i++) {
            y(i, shots_elapsed) = x[det+i] ? 1 : -1;
        }
        shots_elapsed++;
        if (shots_elapsed % 100000 == 0)    std::cout << "shots_elapsed: " << shots_elapsed << "\n";
    };
    generate_syndromes(training_circuit, shots, cb);
    data_matrix.reshape(det, shots_elapsed);
    y.reshape(obs, shots_elapsed);

    ens::Adam opt(1e-3, 512);
    opt.MaxIterations() = shots_elapsed*config.max_epochs;

    model.Train(data_matrix, y, opt, ens::ProgressBar(), ens::PrintLoss());
}

Decoder::result_t
NeuralDecoder::decode_error(const syndrome_t& syndrome) {
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();

    arma::mat data_matrix(det, 1, arma::fill::value(-1.0));
    
    timer.clk_start();

    const std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    for (uint d : detectors)    data_matrix(d, 0) = 1;
    arma::mat pred;
    model.Predict(data_matrix, pred);

    fp_t t = (fp_t)timer.clk_end();

    stim::simd_bits corr(obs);
    corr.clear();
    for (uint i = 0; i < obs; i++) {
        corr[i] ^= (bool)(pred(i, 0) > 0);
    }

    Decoder::result_t res = {
        t,
        corr,
        is_error(corr, syndrome)
    };
    return res;
}
}   // qontra
