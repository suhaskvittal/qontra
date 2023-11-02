/* author: Suhas Vittal
 *  date:   28 August 2023
 * */


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
#endif

using namespace qontra;

stim::Circuit
get_circuit(const schedule_t& sch, fp_t p) {
    const uint n = get_number_of_qubits(sch);

    tables::ErrorAndTiming et;
    et = et * (1000 * p);
    et.e_g1q = 0.0;
    et.e_g2q = 0.0;
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

    uint64_t tshots = 10'000'000;
    uint64_t epochs = 100;

    if (!pp.get_string("asm", asm_file))    return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_float("p", p))  return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_uint64("epochs", epochs);
    pp.get_uint64("tshots", tshots);

    // Get schedule from file.
    schedule_t sch = schedule_from_file(asm_file);
    // Define Decoder.
    stim::Circuit error_model = get_circuit(sch, p);
#ifdef USE_NEURAL_NET
    using namespace mlpack;
    NeuralDecoder dec(error_model);
    dec.model.Add<Linear>(256);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(64);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(error_model.count_observables());
    dec.model.Add<TanH>();
    dec.config.max_epochs = epochs;
    dec.training_circuit = get_circuit(sch, p);
#else
    //MWPMDecoder dec(error_model);
    RestrictionDecoder dec(error_model);
#endif

    // Setup experiment.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
#ifdef DISABLE_MPI
    experiments::G_USE_MPI = false;
#endif
    experiments::memory_params_t params;
    params.shots = shots;
    
#ifdef USE_NEURAL_NET
    std::cout << "starting training...\n";
    dec.train(tshots);
#endif
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
            << res.logical_error_rate;
        for (uint i = 0; i < error_model.count_observables(); i++) {
            out << "," << res.logical_error_rate_by_obs[i];
        }
        out << "," << res.hw_mean << ","
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
