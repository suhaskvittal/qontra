/*
 *  author: Suhas Vittal
 *  date:   8 August 2023
 * */

#ifndef NEURAL_DECODER_h
#define NEURAL_DECODER_h

#define MLPACK_ENABLE_ANN_SERIALIZATION

#include "qontra/decoder.h"

#include <mlpack/methods/ann/ffn.hpp>

#include <vector>

namespace qontra {

class NeuralDecoder : public Decoder {
public:
    NeuralDecoder(DetailedStimCircuit);

    virtual void        train(uint64_t shots, bool verbose=true);
    Decoder::result_t   decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
    virtual void        load_model_from_file(std::string);
    virtual void        save_model_to_file(std::string);

    struct {
        int max_epochs = 100;
    } config;

    mlpack::FFN<mlpack::MeanSquaredError> model;
    DetailedStimCircuit training_circuit;
};

class FragmentedNeuralDecoder : public NeuralDecoder {
public:
    FragmentedNeuralDecoder(DetailedStimCircuit);

    void                train(uint64_t shots, bool verbose=true) override;
    Decoder::result_t   decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
    void                load_model_from_folder(std::string);
    void                save_model_to_folder(std::string);

    // These two functions just call load_model_from_folder and save_model_to_folder.
    void                load_model_from_file(std::string) override;
    void                save_model_to_file(std::string) override;
private:
    std::string get_model_filename(size_t);
    void        set_all_configs(void);

    // The number of backing decoders is equal to the number of observables in the circuit.
    std::vector<uptr<NeuralDecoder>> backing_decoders;
};

std::vector<DetailedStimCircuit> isolate_observables(DetailedStimCircuit);

}   // qontra

#include "neural.inl"

#endif  // NEURAL_DECODER_h
