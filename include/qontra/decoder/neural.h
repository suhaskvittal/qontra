/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#ifndef NEURAL_DECODER_h
#define NEURAL_DECODER_h

#define MLPACK_ENABLE_ANN_SERIALIZATION

#include "qontra/decoder/decoder.h"

#include <mlpack/methods/ann/ffn.hpp>

namespace qontra {

class NeuralDecoder : public Decoder {
public:
    NeuralDecoder(DetailedStimCircuit circ)
        :Decoder(circ, graph::DecodingGraph::Mode::DO_NOT_BUILD),
        training_circuit(circ)
    {}

    void                train(uint64_t shots, bool verbose=true);
    Decoder::result_t   decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
    void                load_model_from_file(std::string);
    void                save_model_to_file(std::string);

    std::string name() override { return "NeuralDecoder"; }

    struct {
        int max_epochs = 100;
    } config;

    mlpack::FFN<mlpack::MeanSquaredError> model;
    DetailedStimCircuit training_circuit;
};

}   // qontra

#endif  // NEURAL_DECODER_h
