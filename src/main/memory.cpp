/* author: Suhas Vittal
 *  date:   28 August 2023
 * */

//#define ARMA_OPENMP_THREADS 32
//#define DISABLE_MPI
//#define USE_NEURAL_NET

#include <decoder/mwpm.h>
#include <decoder/restriction.h>
#include <experiments.h>
#include <parsing/cmd.h>
#include <instruction.h>
#include <tables.h>

#include <filesystem>
#include <fstream>
#include <iostream>

#ifdef USE_NEURAL_NET
#include <decoder/neural.h>
#include <armadillo>
#endif

using namespace qontra;

stim::Circuit
get_circuit(const schedule_t& sch, fp_t p) {
    const uint n = get_number_of_qubits(sch);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;

//  et.e_m1w0 = 0.0;
//  et.e_m0w1 = 0.0;
//  et.e_g1q = 0.0;
    et.e_g2q = 0.0;
    et.e_idle = 0.0;
    et = et * (1000 * p);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    stim::Circuit circ = schedule_to_stim(sch, errors, timing, p);
    return circ;
}

int main(int argc, char* argv[]) {
    int world_rank = 0, world_size = 1;
#ifndef DISABLE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    CmdParser pp(argc, argv);

    std::string asm_file;
    std::string output_file;
    fp_t p;
    uint64_t shots;

    int32_t seed = 0;

    uint64_t tshots = 0;
    uint64_t epochs = 100;
    std::string model_file = "model.bin";

    if (!pp.get_string("asm", asm_file))    return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_float("p", p))  return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_int32("seed", seed);

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
    schedule_t sch = schedule_from_file(asm_file);
    // Define Decoder.
    stim::Circuit error_model = get_circuit(sch, p);
#ifdef USE_NEURAL_NET
    using namespace mlpack;
    NeuralDecoder dec(error_model);
    // Check if model file exists. If so, load it in. 
    // If not, then make and train it.
    std::filesystem::path model_file_path(model_file);
    if (std::filesystem::exists(model_file_path)) {
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
        dec.training_circuit = get_circuit(sch, p);

        std::cout << "starting training...\n";
        dec.train(tshots);

        dec.save_model_to_file(model_file);
    }
#else
    //MWPMDecoder dec(error_model);
    uint64_t detectors_per_round;
    if (!pp.get_uint64("dpr", detectors_per_round)) return 1;
    RestrictionDecoder dec(error_model, detectors_per_round);
#endif
    experiments::memory_params_t params;
    params.shots = shots;
    // Run experiment.
    experiments::memory_result_t res = memory_experiment(&dec, params);

    // Write results to file.
    std::filesystem::path output_path(output_file);
    if (world_rank == 0) {
        std::filesystem::path output_folder(output_path.parent_path());
        safe_create_directory(output_folder);
    }
    bool write_header = !std::filesystem::exists(output_path);
#ifndef DISABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    std::ofstream out(output_path, std::ios::app);
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
        out << std::filesystem::path(asm_file).filename() << ","
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
