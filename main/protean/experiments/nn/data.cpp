/*
 *  author: Suhas Vittal
 *  date:   25 February 2024
 * */

#include <qontra/decoder/mwpm.h>
#include <qontra/decoder/restriction.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#include <mlpack.hpp>

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    vtils::CmdParser pp(argc, argv, 1);
    std::string protean_folder(argv[1]);
    // Parameters for building the training data:
    std::string decoder_name;
    uint64_t    shots = 1'000'000;
    fp_t        p = 1e-3;

    pp.get("decoder", decoder_name, true);
    pp.get("s", shots);
    pp.get("p", p);

    // Make the trace output directory.
    std::string qes_file = protean_folder + "/ext.qes";
    std::string training_data_folder = protean_folder + "/nn/data";
    vtils::safe_create_directory(training_data_folder);

    DetailedStimCircuit circuit = make_circuit(qes_file, p);
    uptr<Decoder> dec = nullptr;
    if (decoder_name == "mwpm") {
        dec = std::make_unique<MWPMDecoder>(circuit);
    } else if (decoder_name == "restriction") {
        dec = std::make_unique<RestrictionDecoder>(circuit);
    } else {
        std::cerr << "Unsupported decoder type: " << decoder_name << std::endl;
    }

    G_FILTER_OUT_SYNDROMES = true;
    G_FILTERING_HAMMING_WEIGHT = 0;
    // Make training data.
    const uint64_t local_shots = shots / world_size + (world_rank==0)*(shots % world_size);
    uint64_t s = 0;

    // Data is column major, so shots is the second index.
    const size_t n_det = circuit.count_detectors();
    const size_t n_obs = circuit.count_observables();
    const size_t n_flags = circuit.flag_detectors.size();

    arma::mat data_matrix(n_flags+n_obs, local_shots, arma::fill::value(-1.0)),
              labels(n_obs, local_shots, arma::fill::zeros);
    generate_syndromes(circuit, shots,
        [&] (shot_payload_t payload)
        {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = payload.syndrome,
                                                    obs = payload.observables;
            // Count number of detectors and flags.
            size_t nf = 0;
            for (uint64_t f : circuit.flag_detectors) {
                nf += syndrome[f];
            }
            if (nf == 0) return;
            size_t hw = syndrome.popcnt();
            if (G_FILTER_OUT_SYNDROMES && hw - nf <= G_FILTERING_HAMMING_WEIGHT) {
                return;
            }
            // Have base decoder decode the syndrome.
            auto res = dec->decode_error(syndrome);
            // Write flags to matrix.
            size_t k = 0;
            for (uint64_t f : circuit.flag_detectors) {
                data_matrix(k, s) = syndrome[f] ? 1 : -1;
                k++;
            }
            for (k = 0; k < n_obs; k++) {
                data_matrix(n_flags+k, s) = res.corr[k] ? 1 : -1;
                labels(k, s) = (obs[k] ^ res.corr[k]) ? 1 : -1;
            }
            s++;
        });
    data_matrix.reshape(n_flags+n_obs, s);
    labels.reshape(n_obs, s);
    // Write the matrices to a file.
    std::string data_file = training_data_folder + "/data." + std::to_string(world_rank) + ".bin";
    std::string labels_file = training_data_folder + "/labels." + std::to_string(world_rank) + ".bin";
    mlpack::data::Save(data_file, data_matrix, true);
    mlpack::data::Save(labels_file, labels, true);

    MPI_Finalize();
    return 0;
}
