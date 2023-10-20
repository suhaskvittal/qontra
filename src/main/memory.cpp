/*
 *  author: Suhas Vittal
 *  date:   28 August 2023
 * */

#define ARMA_OPENMP_THREADS 32

#include "decoder/mwpm.h"
#include "decoder/neural.h"
#include "decoder/restriction.h"
#include "experiments.h"
#include "parsing/cmd.h"
#include "instruction.h"
#include "tables.h"

#include <filesystem>
#include <fstream>
#include <iostream>

#include <omp.h>

#include <armadillo>

using namespace qontra;

#define DISABLE_MPI

stim::Circuit
get_circuit(const schedule_t& sch, fp_t p) {
    const uint n = get_number_of_qubits(sch);

    tables::ErrorAndTiming et;
    et = et * (1000 * p);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    stim::Circuit circ = schedule_to_stim(sch, errors, timing);
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

    uint64_t tshots = 10'000'000;
    uint64_t epochs = 100;

    if (!pp.get_string("asm", asm_file))    return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_float("p", p))  return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_uint64("epochs", epochs);
    pp.get_uint64("tshots", tshots);

    std::cout << "arma config: " << arma::arma_config::mp_threads << "\n";
    std::cout << "arma threads: " << arma::mp_thread_limit::get() << "," << arma::mp_thread_limit::in_parallel() << "\n";
    std::cout << "omp threads: " << omp_get_max_threads() << "\n";

    // Get schedule from file.
    schedule_t sch = schedule_from_file(asm_file);
    // Define Decoder.
    using namespace mlpack;
    stim::Circuit error_model = get_circuit(sch, p);
    NeuralDecoder dec(error_model);
    dec.model.Add<Linear>(256);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(64);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(1);
    dec.model.Add<TanH>();
    dec.config.max_epochs = epochs;
    dec.training_circuit = get_circuit(sch, p);
    /*
    MWPMDecoder dec(error_model);
    RestrictionDecoder dec(error_model);
    */

    // Setup experiment.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
#ifdef DISABLE_MPI
    experiments::G_USE_MPI = false;
#endif
    experiments::memory_params_t params;
    params.shots = shots;
    
    dec.train(tshots);
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
