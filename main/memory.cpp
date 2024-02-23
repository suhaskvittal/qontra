/*
 *  author: Suhas Vittal
 *  date:   21 January 2024
 * */

#include <qontra/decoder/chromobius.h>
#include <qontra/decoder/mwpm.h>
#include <qontra/decoder/neural_assisted.h>
#include <qontra/decoder/pymatching.h>
#include <qontra/decoder/restriction.h>
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

    CmdParser pp(argc, argv);
    std::string HELP = 
        "usage: ./memory --qes <file> --out <file>\n"
        "optional:\n"
        "\t--s <shots, default=1e6>\n"
        "\t--pmin <error-rate, default=1e-3>\n"
        "\t--pmax <error-rate, default=1e-3>\n"
        "\t--steps <# of p, default=1>\n";

    std::string qes_file;
    std::string output_file;

    uint64_t    shots = 1'000'000;
    fp_t        pmin = 1e-3,
                pmax = 1e-3;
    uint64_t    steps = 1;

    pp.get("qes", qes_file, true);
    pp.get("out", output_file, true);

    pp.get("s", shots);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("steps", steps);

    qes::Program<> program = qes::from_file(qes_file);

    DetailedStimCircuit _circuit = make_circuit(qes_file, pmax, true);

    RestrictionDecoder base(_circuit);
    NeuralAssistedDecoder dec(_circuit, &base);

    using namespace mlpack;
    dec.model.Add<Linear>(256);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(64);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(_circuit.count_observables());
    dec.model.Add<TanH>();

    dec.config.max_epochs = 1000;

    dec.train(1'000'000);

    fp_t p = pmin;
    while (p <= 1.1*pmax) {
        // Load model from file and run memory experiment.
        DetailedStimCircuit circuit = make_circuit(qes_file, p, true);
        dec.set_circuit(circuit);

        memory_config_t config;
        config.shots = shots;
        auto res = memory_experiment(&dec, config);

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
                << res.logical_error_rate;
            for (fp_t x : res.logical_error_rate_by_obs) {
                fout << "," << x;
            }
            fout << std::endl;
        }
        p = pow(M_E, log(p) + (log(pmax)-log(pmin))/(steps-1));
    }
    MPI_Finalize();
    return 0;
}
