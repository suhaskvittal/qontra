/*
 *  author: Suhas Vittal
 *  date:   21 January 2024
 * */

#include <qontra/decoder/pymatching.h>
#include <qontra/decoder/restriction.h>
#include <qontra/ext/stim.h>

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
    configure_optimal_batch_size();
    // Initialize MPI
    MPI_Init(NULL, NULL);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser pp(argc, argv, 1);

    std::string qes_file(argv[1]);
    std::string output_file(argv[2]);
    std::string decoder;

    if (world_rank == 0) {
        std::cout << "reading " << qes_file << ", writing to " << output_file << std::endl;
    }

    uint64_t    errors_until_stop = 40;
    fp_t        pmin = 5e-4,
                pmax = 3e-3;
    uint64_t    step_size = 1;

    pp.get("decoder", decoder, true);

    pp.get("e", errors_until_stop);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("step-size", step_size);

    G_BASE_SEED = 1000;

    // Initialize error and timing tables.
    qes::Program<> program = qes::from_file(qes_file);

    DetailedStimCircuit base_circuit = make_default_circuit(program, pmax);

    uptr<Decoder> dec;
    if (decoder == "mwpm") {
        dec = std::make_unique<PyMatching>(base_circuit);
    } else {
        dec = std::make_unique<RestrictionDecoder>(base_circuit);
    }

    memory_config_t config;
    config.errors_until_stop = errors_until_stop;
    fp_t p = pmin;
    while (p <= 1.1*pmax) {
        // Load model from file and run memory experiment.
        DetailedStimCircuit circuit = make_default_circuit(program, p);
        dec->set_circuit(circuit);
        auto res = run_memory_with_generated_syndromes(dec.get(), config);
        // Write result to file.
        if (world_rank == 0) {
            std::cout << "[ pr_base_memory ] Writing p = " << p << std::endl;
            bool write_header = !file_exists(output_file);
            if (!file_exists(get_parent_directory(output_file.c_str()))) {
                safe_create_directory(get_parent_directory(output_file.c_str()));
            }
            std::ofstream fout(output_file, std::ios::app);
            if (write_header) {
                fout << "physical error rate,logical error rate" << std::endl;
            }
            fout << p << ","
                << res.logical_error_rate;
            fout << std::endl;
        }
        fp_t gran = floor(log10(p));
        p += step_size*pow(10, gran);
    }
    MPI_Finalize();
    return 0;
}

