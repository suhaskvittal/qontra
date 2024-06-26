/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include <qontra/decoder/mwpm.h>
#include <qontra/decoder/neural_assisted.h>
#include <qontra/decoder/restriction.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#include <fstream>
#include <iostream>

#include <mpi.h>

using namespace qontra;
using namespace vtils;

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(NULL, NULL);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser pp(argc, argv, 1);

    std::string protean_folder(argv[1]);
    std::string qes_file = protean_folder + "/ext.qes";
    std::string output_file = protean_folder + "/output/nn_memory.csv";

    std::string decoder_name;

    uint64_t    shots = 1'000'000;
    fp_t        pmin = 1e-4,
                pmax = 1e-3;
    size_t      steps = 10;

    pp.get("decoder", decoder_name, true);

    pp.get("s", shots);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("steps", steps);

    G_BASE_SEED = 1000;

    qes::Program<> program = qes::from_file(qes_file);
    DetailedStimCircuit _circuit = make_circuit(program, pmax);

    uptr<Decoder> base = nullptr;
    if (decoder_name == "mwpm") {
        base = std::make_unique<MWPMDecoder>(_circuit);
    } else if (decoder_name == "restriction") {
        base = std::make_unique<RestrictionDecoder>(_circuit);
    } else {
        std::cerr << "Unsupported decoder type: " << decoder_name << std::endl;
    }
    uptr<NeuralAssistedDecoder> dec = std::make_unique<NeuralAssistedDecoder>(_circuit, base.get());
    // Load model for neural assisted decoder.
    dec->load_model_from_file(protean_folder + "/nn/model.bin");

    fp_t p = pmin;
    while (p <= 1.01*pmax) {
        // Load model from file and run memory experiment.
        DetailedStimCircuit circuit = make_circuit(program, p);
        dec->set_circuit(circuit);

        memory_config_t config;
        config.shots = shots;
        auto res = memory_experiment(dec.get(), config);

        // Write result to file.
        if (world_rank == 0) {
            bool write_header = !file_exists(output_file);
            if (!file_exists(get_parent_directory(output_file.c_str()))) {
                safe_create_directory(get_parent_directory(output_file.c_str()));
            }
            std::ofstream fout(output_file, std::ios::app);
            if (write_header) {
                fout << "physical error rate,shots,logical error rate" << std::endl;
            }
            fout << p << ","
                << shots << ","
                << res.logical_error_rate;
            fout << std::endl;
        }
        p = pow(M_E, log(p) + (log(pmax)-log(pmin))/(steps-1));
    }
    MPI_Finalize();
    return 0;
}
