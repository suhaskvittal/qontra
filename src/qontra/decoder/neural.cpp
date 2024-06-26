/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#include "qontra/decoder/neural.h"
#include "qontra/experiments.h"

#include <vtils/filesystem.h>

namespace qontra {

static const uint64_t PER_HW_LIMIT = 100'000;

NeuralDecoder::NeuralDecoder(const DetailedStimCircuit& circuit)
    :Decoder(circuit),
    model(),
    training_circuit(circuit)
{}

void
NeuralDecoder::train(uint64_t shots, bool verbose) {
    const size_t n_det = circuit.count_detectors();
    const size_t n_obs = circuit.count_observables();

    arma::mat data_matrix(n_det, shots, arma::fill::value(-1.0));
    arma::mat y(n_obs, shots, arma::fill::zeros);

    uint64_t shots_elapsed = 0;
    // We will train the model every batch.
    std::map<size_t, uint64_t> hw_freq;
    generate_syndromes(training_circuit, shots,
        [&] (shot_payload_t payload)
        {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = payload.syndrome;
            stim::simd_bits_range_ref<SIMD_WIDTH> obs = payload.observables;

            const size_t hw = syndrome.popcnt();
            if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
                return;
            }
            if (hw_freq[hw]++ >= PER_HW_LIMIT)  return;

            std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
            for (uint64_t d : detectors) {
                data_matrix(d, shots_elapsed) = 1;
            }
            for (size_t i = 0; i < n_obs; i++) {
                y(i, shots_elapsed) = obs[i] ? 1 : -1;
            }
            shots_elapsed++;

            if (shots_elapsed % 100'000 == 0 && verbose) {
                std::cout << "[ NeuralDecoder ] shots_elapsed: " << shots_elapsed << std::endl;
            }
        });
    data_matrix.reshape(n_det, shots_elapsed);
    y.reshape(n_obs, shots_elapsed);

    ens::Adam opt(1e-3, 512);
    opt.MaxIterations() = shots_elapsed*config.max_epochs;

    model.Train(data_matrix, y, opt, ens::ProgressBar(), ens::PrintLoss());
}

Decoder::result_t
NeuralDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_det = circuit.count_detectors();
    const size_t n_obs = circuit.count_observables();

    arma::mat data_matrix(n_det, 1, arma::fill::value(-1.0));
    
    timer.clk_start();

    const std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
    for (uint64_t d : detectors) data_matrix(d, 0) = 1;
    arma::mat pred;
    model.Predict(data_matrix, pred);

    fp_t t = static_cast<fp_t>(timer.clk_end());

    stim::simd_bits<SIMD_WIDTH> corr(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        corr[i] ^= (bool)(pred(i, 0) > 0);
    }

    Decoder::result_t res = {
        t,
        corr
    };
    return res;
}

void
NeuralDecoder::save_model_to_file(std::string fname) {
    mlpack::data::Save(fname, "model", model, false);
}

void
NeuralDecoder::load_model_from_file(std::string fname) {
    mlpack::data::Load(fname, "model", model);
}

}   // qontra
