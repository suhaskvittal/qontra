/*
 *  author: Suhas Vittal
 *  date:   28 August 2023
 * */

#include "decoder/mwpm.h"
#include "decoder/pymatching.h"
#include "mldh/astrea.h"
#include "mldh/block.h"
#include "experiments.h"
#include "parsing/cmd.h"
#include "instruction.h"
#include "tables.h"

#include <filesystem>
#include <fstream>
#include <iostream>

#include <math.h>
#include <mpi.h>

using namespace qontra;
using namespace mldh;

#define BMEM_SQR(x)          ((x) * (x))
#define BMEM_MEAN(x, n)      ((fp_t)(x)) / ((fp_t)(n))
#define BMEM_STD(x, x2, n)   sqrt(BMEM_MEAN(x2, n) - BMEM_SQR(BMEM_MEAN(x, n)))

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
    uint blk;   // block size

    std::string output_file;
    uint64_t shots;

    if (!pp.get_uint32("d", d))             return 1;
    if (!pp.get_float("p", p))              return 1;
    if (!pp.get_uint32("r", r))             return 1;
    if (!pp.get_uint32("blk", blk))         return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_uint64("shots", shots))     return 1;

    // Define circuits.
    stim::CircuitGenParameters circ_params(r, d, "rotated_memory_z");
    set_error_rate_to(p, circ_params);
    stim::Circuit circuit = stim::generate_surface_code_circuit(circ_params).circuit;

    AstreaDecoder<PyMatching, 8> base_dec(circuit);
    BlockDecoder dec(circuit, &base_dec, blk);

    // Setup experiment.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
    experiments::G_FILTERING_HAMMING_WEIGHT = 2;
    experiments::memory_params_t params;
    params.shots = shots;
    
    // Run experiment.
    experiments::memory_params_t caching_params;
    caching_params.shots = world_size * 10;
    memory_experiment(&dec, caching_params);

    dec.reset_stats();
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
    
    // Accumulate custom statistics.
    uint64_t total_number_of_blocks;
    uint64_t total_number_sqr_of_blocks;
    uint64_t max_number_of_blocks;
    uint64_t min_number_of_blocks;
    uint64_t true_shots;

    uint64_t total_hw_in_block;
    uint64_t total_hw_sqr_in_block;
    uint64_t max_hw_in_block;

    uint64_t total_blk_hw_above_th;

    MPI_Reduce(&dec.total_number_of_blocks, &total_number_of_blocks, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_number_sqr_of_blocks, &total_number_sqr_of_blocks, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.max_number_of_blocks, &max_number_of_blocks, 1,
            MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.min_number_of_blocks, &min_number_of_blocks, 1,
            MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_shots_evaluated, &true_shots, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&dec.total_hw_in_block, &total_hw_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_hw_sqr_in_block, &total_hw_sqr_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.max_hw_in_block, &max_hw_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&dec.total_blk_hw_above_th, &total_blk_hw_above_th, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        fp_t mean_blocks = BMEM_MEAN(total_number_of_blocks, true_shots);
        fp_t std_blocks = BMEM_STD(total_number_of_blocks,
                                total_number_sqr_of_blocks,
                                true_shots);
        fp_t mean_blkhw = BMEM_MEAN(total_hw_in_block, total_number_of_blocks);
        fp_t std_blkhw = BMEM_STD(total_hw_in_block,
                                    total_hw_sqr_in_block,
                                    total_number_of_blocks);
        fp_t prob_blkhw_gtth = BMEM_MEAN(total_blk_hw_above_th, total_number_of_blocks);

        if (write_header) {
            // Write the header.
            out << "Distance,"
                    << "Rounds,"
                    << "Physical Error Rate,"
                    << "Block Width,"
                    << "Shots,"
                    << "Logical Error Probability,"
                    << "Hamming Weight Mean,"
                    << "Hamming Weight Std,"
                    << "Hamming Weight Max,"
                    << "Time Mean,"
                    << "Time Std,"
                    << "Time Max,"
                    << "Block Mean,"
                    << "Block Std,"
                    << "Block Min,"
                    << "Block Max,"
                    << "Block Hamming Weight Mean,"
                    << "Block Hamming Weight Std,"
                    << "Block Hamming Weight Max,"
                    << "Prob Block Hamming Weight GT " << dec.config.blocking_threshold << "\n";
        }
        out << d << ","
            << r << ","
            << p << ","
            << blk << ","
            << shots << ","
            << res.logical_error_rate << ","
            << res.hw_mean << ","
            << res.hw_std << ","
            << res.hw_max << ","
            << res.t_mean << ","
            << res.t_std << ","
            << res.t_max << ","
            << mean_blocks << ","
            << std_blocks << ","
            << min_number_of_blocks << ","
            << max_number_of_blocks << ","
            << mean_blkhw << ","
            << std_blkhw << ","
            << max_hw_in_block << ","
            << prob_blkhw_gtth << "\n";
    }
    MPI_Finalize();
}
