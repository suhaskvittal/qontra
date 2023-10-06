/*
 *  author: Suhas Vittal
 *  date:   28 August 2023
 * */

#include "decoder/mwpm.h"
#include "decoder/neural.h"
#include "experiments.h"
#include "parsing/cmd.h"
#include "instruction.h"
#include "tables.h"

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CmdParser pp(argc, argv);

    std::string asm_file;
    std::string output_file;
    fp_t p;
    uint64_t shots;

    if (!pp.get_string("asm", asm_file))    return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_float("p", p))  return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    // Get schedule from file.
    schedule_t sch = schedule_from_file(asm_file);
    const uint n = get_number_of_qubits(sch);

    // Define error model.
    tables::ErrorAndTiming et;
    et = et * (p * 1000);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);

    // Define Decoder.
    using namespace mlpack;
    stim::Circuit error_model = schedule_to_stim(sch, errors, timing);
    NeuralDecoder dec(error_model);
    dec.model.Add<Linear>(256);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(64);
    dec.model.Add<TanH>();
    dec.model.Add<Linear>(1);
    dec.model.Add<TanH>();
    dec.config.max_epochs = 25;
    dec.train(10*shots);
    /*
    MWPMDecoder dec(error_model);
    */

    // Setup experiment.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
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
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Finalize();
}
