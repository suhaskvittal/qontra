/*
 *  author: Suhas Vittal
 *  date:   28 August 2023
 * */

#include "decoder/mwpm.h"
#include "decoder/window.h"
#include "experiments.h"
#include "parsing/cmd.h"
#include "instruction.h"
#include "tables.h"

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace qontra;

void set_error_rate_to(fp_t p, stim::CircuitGenParameters& params) {
    params.before_measure_flip_probability = p;
    params.after_reset_flip_probability = p;
    params.after_clifford_depolarization = p;
    params.before_round_data_depolarization = p;
}

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CmdParser pp(argc, argv);

    uint d; // Distance
    fp_t p; // Error rate
    uint r; // Rounds
    uint w; // Window size
    uint c; // Commit window

    std::string output_file;
    uint64_t shots;

    if (!pp.get_uint32("d", d))             return 1;
    if (!pp.get_float("p", p))              return 1;
    if (!pp.get_uint32("r", r))             return 1;
    if (!pp.get_uint32("w", w))             return 1;
    if (!pp.get_uint32("c", c))             return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_uint64("shots", shots))     return 1;

    // Define circuits.
    stim::CircuitGenParameters circ_params(r, d, "rotated_memory_z");
    set_error_rate_to(p, circ_params);
    stim::Circuit circuit = stim::generate_surface_code_circuit(circ_params).circuit;
    
    stim::CircuitGenParameters window_circ_params(w+1, d, "rotated_memory_z");
    set_error_rate_to(p, window_circ_params);
    stim::Circuit window_circuit = stim::generate_surface_code_circuit(window_circ_params).circuit;

    // Define Decoder.
    const uint dpr = (d*d - 1) >> 1;

    MWPMDecoder base_dec(window_circuit);
    WindowDecoder dec(circuit, &base_dec, c, dpr);

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
            out << "Distance,"
                    << "Rounds,"
                    << "Physical Error Rate,"
                    << "Window Read:Commit,"
                    << "Shots,"
                    << "Logical Error Probability,"
                    << "Hamming Weight Mean,"
                    << "Hamming Weight Std,"
                    << "Hamming Weight Max,"
                    << "Time Mean,"
                    << "Time Std,"
                    << "Time Max\n";
        }
        out << d << ","
            << r << ","
            << p << ","
            << w << ":" << c << ","
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
