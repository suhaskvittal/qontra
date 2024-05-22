/*
 *  author: Suhas Vittal
 *  date:   21 January 2024
 * */

#include <qontra/decoder/chromobius.h>
#include <qontra/decoder/mwpm.h>
#include <qontra/decoder/pymatching.h>
#include <qontra/decoder/restriction.h>
#include <qontra/decoder/mobius.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>

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

    CmdParser pp(argc, argv, 2);
    std::string HELP = 
        "usage: ./memory <qes-file> <output-file>\n"
        "optional:\n"
        "\t--e <errors until stopping, default=40>\n"
        "\t--pmin <error-rate, default=1e-3>\n"
        "\t--pmax <error-rate, default=1e-3>\n"
        "\t--step-size <default=1>\n";

    std::string qes_file(argv[1]);
    std::string output_file(argv[2]);

    uint64_t    errors_until_stop = 40;
    fp_t        pmin = 1e-3,
                pmax = 1e-3;
    uint64_t    step_size = 1;

    pp.get("e", errors_until_stop);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("step-size", step_size);

    qes::Program<> program = qes::from_file(qes_file);

    DetailedStimCircuit _circuit = make_default_circuit(qes_file, pmax, true, "cap");
    RestrictionDecoder dec(_circuit);

    memory_config_t config;
    config.errors_until_stop = errors_until_stop;

    fp_t p = pmin;
    while (p <= 1.1*pmax) {
        // Load model from file and run memory experiment.
        DetailedStimCircuit circuit = make_default_circuit(qes_file, p, true, "cap");
        dec.set_circuit(circuit);

        auto res = run_memory_with_generated_syndromes(&dec, config);
        // Write result to file.
        if (world_rank == 0) {
            std::cout << "Writing p = " << p << std::endl;
            bool write_header = false;
            if (!file_exists(get_parent_directory(output_file.c_str()))) {
                safe_create_directory(get_parent_directory(output_file.c_str()));
                write_header = true;
            }
            std::ofstream fout(output_file, std::ios::app);
            if (write_header) {
                fout << "qes file,physical error rate,logical error rate" << std::endl;
            }
            fout << get_basename(qes_file) << ","
                << p << ","
                << res.logical_error_rate;
            for (fp_t x : res.logical_error_rate_by_obs) {
                fout << "," << x;
            }
            fout << std::endl;
        }
        fp_t gran = floor(log10(p));
        p += step_size*pow(10, gran);
    }
    MPI_Finalize();
    return 0;
}
