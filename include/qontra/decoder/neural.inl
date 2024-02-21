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

}   // qontra
