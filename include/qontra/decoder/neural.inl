/*
 *  author: Suhas Vittal
 *  date:   15 January 2024
 * */

#include <vtils/filesystem.h>

namespace qontra {

inline void
NeuralDecoder::save_model_to_file(std::string fname) {
    mlpack::data::Save(fname, "model", model, false);
}

inline void
NeuralDecoder::load_model_from_file(std::string fname) {
    mlpack::data::Load(fname, "model", model);
}

inline void
FragmentedNeuralDecoder::train(uint64_t shots, bool verbose) {
    set_all_configs();
    for (size_t i = 0; i < backing_decoders.size(); i++) {
        if (verbose) {
            std::cout << "[ FragmentedNeuralDecoder, NN " << i << " ]" << std::endl;
        }
        backing_decoders[i]->train(shots, verbose);
    }
}

inline Decoder::result_t
FragmentedNeuralDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    stim::simd_bits<SIMD_WIDTH> corr(backing_decoders.size());
    fp_t t = 0.0;
    for (size_t i = 0; i < backing_decoders.size(); i++) {
        Decoder::result_t res = backing_decoders[i]->decode_error(syndrome);
        // Now, we need to update corr.
        t += res.exec_time;
        corr[i] = res.corr[0];
    }
    return { t, corr };
}

inline void
FragmentedNeuralDecoder::load_model_from_folder(std::string folder) {
    for (size_t i = 0; i < backing_decoders.size(); i++) {
        backing_decoders[i]->load_model_from_file(folder + "/" + get_model_filename(i));
    }
}

inline void
FragmentedNeuralDecoder::save_model_to_folder(std::string folder) {
    vtils::safe_create_directory(folder);
    for (size_t i = 0; i < backing_decoders.size(); i++) {
        backing_decoders[i]->save_model_to_file(folder + "/" + get_model_filename(i));
    }
}

inline void
FragmentedNeuralDecoder::load_model_from_file(std::string folder) {
    load_model_from_folder(folder);
}

inline void
FragmentedNeuralDecoder::save_model_to_file(std::string folder) {
    save_model_to_folder(folder);
}

inline std::string
FragmentedNeuralDecoder::get_model_filename(size_t k) {
    return "model_" + std::to_string(k) + ".bin";
}

inline void
FragmentedNeuralDecoder::set_all_configs() {
    // Copy data from this object to all the backing decoders.
    std::vector<DetailedStimCircuit> training_circuits = isolate_observables(training_circuit);
    for (size_t i = 0; i < backing_decoders.size(); i++) {
        backing_decoders[i]->config = config;
        backing_decoders[i]->training_circuit = training_circuits[i];
        backing_decoders[i]->model = model;
    }
}

}   // qontra
