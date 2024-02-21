/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#include "qontra/decoder/neural_assisted.h"
#include "qontra/experiments.h"

namespace qontra {

NeuralAssistedDecoder::NeuralAssistedDecoder(const DetailedStimCircuit& circuit, Decoder* d_p)
    :NeuralDecoder(circuit),
    dec_p(d_p)
{}

void
NeuralAssistedDecoder::train(uint64_t shots, bool verbose) {
    const size_t n_obs = circuit.count_observables();
    const size_t n_flags = circuit.flag_detectors.size();
    const size_t n = n_flags + n_obs;

    arma::mat data_matrix(n, shots, arma::fill::value(-1.0));
    arma::mat y(n_obs, shots, arma::fill::zeros);

    uint64_t shots_elapsed = 0;
    generate_syndromes(training_circuit, shots,
        [&] (shot_payload_t payload)
        {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = payload.syndrome;
            stim::simd_bits_range_ref<SIMD_WIDTH> obs = payload.observables;
            // Split the detectors into parity and flags.
            std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
            std::vector<uint64_t> flags = get_flags(detectors);
            if (detectors.empty() || flags.empty()) return;
            auto res = dec_p->decode_error(syndrome);
            for (uint64_t f : flags) {
                data_matrix(f, shots_elapsed) = 1;
            }
            for (size_t i = 0; i < n_obs; i++) {
                data_matrix(n_flags+i, shots_elapsed) = res.corr[i] ? 1 : -1;
                y(i, shots_elapsed) = (res.corr[i] ^ obs[i]) ? 1 : -1;
            }
            shots_elapsed++;

            if (shots_elapsed % 100'000 == 0 && verbose) {
                std::cout << "[ NeuralAssistedDecoder ] shots_elapsed: " << shots_elapsed << std::endl;
            }
        });
    data_matrix.reshape(n, shots_elapsed);
    y.reshape(n_obs, shots_elapsed);

    std::cout << "training data size: " << shots_elapsed << "\n";

    ens::Adam opt(1e-3, 512);
    opt.MaxIterations() = shots_elapsed*config.max_epochs;

    model.Train(data_matrix, y, opt, ens::ProgressBar(), ens::PrintLoss());
}

Decoder::result_t
NeuralAssistedDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();
    const size_t n_flags = circuit.flag_detectors.size();
    const size_t n = n_flags + n_obs;

    arma::mat data_matrix(n, 1, arma::fill::value(-1.0));
    
    timer.clk_start();

    std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
    std::vector<uint64_t> flags = get_flags(detectors);
    // Get initial result from the decode.
    auto _res = dec_p->decode_error(syndrome);
    if (detectors.empty() || flags.empty()) return _res;

    for (uint64_t f : flags) {
        data_matrix(f, 0) = 1;
    }
    for (size_t i = 0; i < n_obs; i++) {
        data_matrix(n_flags+i, 0) = _res.corr[i] ? 1 : -1;
    }

    arma::mat pred;
    model.Predict(data_matrix, pred);

    fp_t t = static_cast<fp_t>(timer.clk_end());

    stim::simd_bits<SIMD_WIDTH> corr(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        corr[i] ^= (bool)(pred(i, 0) > 0);
    }

    Decoder::result_t res = {
        t,
        corr ^ _res.corr
    };
    return res;
}

}   // qontra
