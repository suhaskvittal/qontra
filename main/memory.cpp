/* 
 *  author: Suhas Vittal
 *  date:   28 August 2023
 * */

//#define ARMA_OPENMP_THREADS 32
//#define DISABLE_MPI
//#define USE_NEURAL_NET

#include <qontra/decoder/mwpm.h>

#ifdef QONTRA_PYMATCHING_ENABLED
#include <qontra/decoder/pymatching.h>
#endif

#ifdef QONTRA_CHROMOBIUS_ENABLED
#include <qontra/decoder/chromobius.h>
#endif

#include <qontra/experiments.h>
#include <qontra/ext/stim.h>
#include <qontra/ext/qes.h>
#include <qontra/tables.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#include <fstream>
#include <iostream>

#ifdef USE_NEURAL_NET
#include <decoder/neural.h>
#include <armadillo>
#endif

#include <mpi.h>

using namespace qontra;

enum class model_style { cc, pheno, circuit };

DetailedStimCircuit
get_circuit(const qes::Program<>& program, fp_t p, model_style m=model_style::circuit) {
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;

    // TEMPORARY MODIFICATIONS START
    et.e_idle = 0.0;
    // TEMPORARY MODIFICATIONS END

    if (m == model_style::pheno) {
        et.e_g1q = 0.0;
        et.e_g2q = 0.0;
    }
    if (m == model_style::cc) {
        et.e_g1q = 0.0;
        et.e_g2q = 0.0;
        et.e_m1w0 = 0.0;
        et.e_m0w1 = 0.0;
    }
    et = et * (1000 * p);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing);
}

int main(int argc, char* argv[]) {
    int world_rank = 0, world_size = 1;
#ifndef DISABLE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    vtils::CmdParser pp(argc, argv);

    std::string qes_file;
    std::string output_file;
    fp_t p;
    uint64_t shots;

    int32_t seed = 0;
    model_style m = model_style::circuit;

    uint64_t tshots = 0;
    uint64_t epochs = 100;
    std::string model_file = "model.bin";

    DetailedStimCircuit error_model;
    if (pp.get_string("stim", qes_file)) {
        FILE* fin = fopen(qes_file.c_str(), "r");
        error_model = stim::Circuit::from_file(fin);
        fclose(fin);
    } else if (!pp.get_string("qes", qes_file)) {
        return 1;
    }

    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_float("p", p))  return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_int32("seed", seed);

    if (pp.option_set("pheno")) m = model_style::pheno;
    if (pp.option_set("cc")) m = model_style::cc;

    experiments::G_BASE_SEED = seed;

    pp.get_uint64("epochs", epochs);
    pp.get_uint64("tshots", tshots);
    pp.get_string("model-file", model_file);

    // Setup experiment settings.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
#ifdef DISABLE_MPI
    experiments::G_USE_MPI = false;
#endif
    // Get schedule from file.
    if (!pp.option_set("stim")) {
        qes::Program<> program = qes::from_file(qes_file);
        error_model = get_circuit(program, p, m);
    }
#ifdef USE_NEURAL_NET
    using namespace mlpack;
    NeuralDecoder dec(error_model);
    // Check if model file exists. If so, load it in. 
    // If not, then make and train it.
    if (file_exists(model_file)) {
        dec.load_model_from_file(model_file);
    } else {
        dec.model.Add<Linear>(256);
        dec.model.Add<TanH>();
        dec.model.Add<Linear>(64);
        dec.model.Add<TanH>();
        dec.model.Add<Linear>(error_model.count_observables());
        dec.model.Add<TanH>();
    }

    if (tshots > 0) {
        dec.config.max_epochs = epochs;
        dec.training_circuit = get_circuit(sch, p, m);

        std::cout << "starting training...\n";
        dec.train(tshots);

        dec.save_model_to_file(model_file);
    }
#else

#ifdef QONTRA_CHROMOBIUS_ENABLED
//  PyMatching dec(error_model);
    Chromobius dec(error_model);
#else
    MWPMDecoder dec(error_model);
#endif

#endif

    experiments::memory_params_t params;
    params.shots = shots;
    // Run experiment.
    experiments::memory_result_t res = memory_experiment(&dec, params);

    // Write results to file.
    if (world_rank == 0) {
        vtils::safe_create_directory(vtils::get_parent_directory(output_file.c_str()));
    }
    bool write_header = !vtils::file_exists(output_file);
#ifndef DISABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    std::ofstream out(output_file, std::ios::app);
    if (world_rank == 0) {
        if (write_header) {
            // Write the header.
            out << "ASM File,"
                    << "Physical Error Rate,"
                    << "Shots,"
                    << "Logical Error Probability,"
                    << "Hamming Weight Mean,"
                    << "Hamming Weight Std,"
                    << "Hamming Weight Max,"
                    << "Time Mean,"
                    << "Time Std,"
                    << "Time Max\n";
        }
        out << qes_file << ","
            << p << ","
            << shots << ","
            << res.logical_error_rate << ","
            << res.hw_mean << ","
            << res.hw_std << ","
            << res.hw_max << ","
            << res.t_mean << ","
            << res.t_std << ","
            << res.t_max << "\n";
    }
#ifndef DISABLE_MPI
    MPI_Finalize();
#endif
}
