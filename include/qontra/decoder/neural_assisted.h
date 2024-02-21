/*
 *  author: Suhas Vittal
 *  date:   19 February 2024
 * */

#ifndef NEURAL_ASSISTED_h
#define NEURAL_ASSISTED_h

#include "qontra/decoder/neural.h"

namespace qontra {

class NeuralAssistedDecoder : public NeuralDecoder {
public:
    NeuralAssistedDecoder(const DetailedStimCircuit&, Decoder*);

    void                train(uint64_t, bool verbose=true) override;
    Decoder::result_t   decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
private:
    std::vector<uint64_t> get_flags(std::vector<uint64_t>&);

    Decoder* dec_p;
};

}   // qontra

#include "neural_assisted.inl"

#endif  // NEURAL_ASSISTED_h
