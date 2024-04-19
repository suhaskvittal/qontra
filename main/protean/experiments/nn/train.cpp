/*
 *  author: Suhas Vittal
 *  date:   26 February 2024
 * */

#define MLPACK_ENABLE_ANN_SERIALIZATION

#include <vtils/filesystem.h>
#include <mlpack.hpp>

#include <string>

#include <stdlib.h>

inline std::string
get_data_file(std::string folder, int k) {
    return folder + "/nn/data/data." + std::to_string(k) + ".bin";
}

inline std::string
get_labels_file(std::string folder, int k) {
    return folder + "/nn/data/labels." + std::to_string(k) + ".bin";
}

using namespace mlpack;

int main(int argc, char* argv[]) {
    std::string protean_folder(argv[1]);

    std::string model_file = protean_folder + "/nn/model.bin";
    int epochs = atoi(argv[2]);

    arma::mat data_matrix, labels;
    data::Load(get_data_file(protean_folder, 0), data_matrix, true);
    data::Load(get_labels_file(protean_folder, 0), labels, true);

    FFN<MeanSquaredError> model;
    // Setup model layers.
    model.Add<Linear>(256);
    model.Add<TanH>();
    model.Add<Linear>(64);
    model.Add<TanH>();
    model.Add<Linear>(labels.n_rows);
    model.Add<TanH>();
    // Now, we need to train the model. First, we should accumulate the data.
    // If there are other data files, concatenate them to the existing matrices.
    int k = 1;
    while (true) {
        if (!vtils::file_exists(get_data_file(protean_folder, k))) break;

        arma::mat _x, _y;
        data::Load(get_data_file(protean_folder, k), _x, true);
        data::Load(get_labels_file(protean_folder, k), _y, true);
        data_matrix = arma::join_rows(data_matrix, _x);
        labels = arma::join_rows(labels, _y);
        k++;
    }
    // Train the model.
    ens::Adam opt(1e-3, 512);
    opt.MaxIterations() = data_matrix.n_cols*epochs;

    model.Train(data_matrix, labels, opt, ens::ProgressBar(), ens::PrintLoss());
    mlpack::data::Save(model_file, "model", model, false);

    return 0;
}
