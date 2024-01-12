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

#include <stdlib.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    vtils::CmdParser pp(argc, argv);
    std::string HELP = 
        "usage: ./pr_nn_train --qes <file> --out <model-file>\n"
        "\toptional: --s <shots, default=1e6> --e <epochs> --p <error-rate, default=1e-3>";

    std::string qes_file;
    std::string model_file;

    uint64_t    shots = 1'000'000;
    uint64_t    epochs = 100;

    fp_t        p = 1e-3;

    pp.get_string("qes", qes_file, true);
    pp.get_string("out", model_file, true);

    pp.get_uint64("s", shots);
    pp.get_uint64("e", epochs);
    pp.get_float("p", p);

    experiments::G_USE_MPI = false;
    experiments::G_BASE_SEED = rand();

    DetailedStimCircuit circuit = make_circuit(qes_file, p);
    NeuralDecoder dec(circuit);
    dec.config.max_epochs = epochs;
    dec.training_circuit = circuit;

    if (vtils::file_exists(model_file)) {
        dec.load_model_from_file(model_file);
    } else {
        using namespace mlpack;
        dec.model.Add<Linear>(256);
        dec.model.Add<TanH>();
        dec.model.Add<Linear>(64);
        dec.model.Add<TanH>();
        dec.model.Add<Linear>(circuit.count_observables());
        dec.model.Add<TanH>();
    }

    dec.train(shots);
    dec.save_model_to_file(model_file);
    return 0;
}
