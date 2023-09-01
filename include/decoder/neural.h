/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#ifndef NEURAL_DECODER_h
#define NEURAL_DECODER_h

#include "decoder/decoder.h"
#include "experiments.h"

#include <mlpack/methods/ann/ffn.hpp>

namespace qontra {

class NeuralDecoder : public Decoder {
public:
    NeuralDecoder(const stim::Circuit& circ)
        :Decoder(circ, graph::DecodingGraph::Mode::LOW_MEMORY)
    {}

    void                train(uint64_t shots, bool verbose=true);
    Decoder::result_t   decode_error(const syndrome_t&) override;

    std::string         name() override { return "NeuralDecoder"; }

    mlpack::FFN<mlpack::MeanSquaredError> model;
};

}   // qontra

#endif  // NEURAL_DECODER_h
