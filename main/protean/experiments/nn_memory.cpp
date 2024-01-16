/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include "memory.h"

#include <qontra/decoder/neural.h>
#include <qontra/experiments.h>
#include <qontra/ext/stim.h>
#include <qontra/tables.h>

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

    CmdParser pp(argc, argv);
    std::string HELP = 
        "usage: ./pr_nn_memory --qes <file> --out <file> --model <file>\n"
        "\toptional: --s <shots, default=1e6> --p <error-rate, default=1e-3>";

    std::string qes_file;
    std::string model_file;
    std::string output_file;

    uint64_t    shots = 1'000'000;
    fp_t        p = 1e-3;

    pp.get_string("qes", qes_file, true);
    pp.get_string("model", model_file, true);
    pp.get_string("out", output_file, true);

    pp.get_uint64("s", shots);
    pp.get_float("p", p);

    experiments::G_BASE_SEED = 1;

    // Load model from file and run memory experiment.
    DetailedStimCircuit circuit = make_circuit(qes_file, p);
    NeuralDecoder dec(circuit);
    
    dec.load_model_from_file(model_file);

    experiments::memory_params_t params;
    params.shots = shots;
    auto res = memory_experiment(&dec, params);

    // Write result to file.
    if (world_rank == 0) {
        bool write_header = false;
        if (!file_exists(get_parent_directory(output_file.c_str()))) {
            safe_create_directory(get_parent_directory(output_file.c_str()));
            write_header = true;
        }
        std::ofstream fout(output_file, std::ios::app);
        if (write_header) {
            fout << "qes file,physical error rate,shots,logical error rate" << std::endl;
        }
        fout << get_basename(qes_file) << ","
            << p << ","
            << shots << ","
            << res.logical_error_rate << std::endl;
    }
    MPI_Finalize();
    return 0;
}
