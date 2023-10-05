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

void
build_folder(std::string folder) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_rank == 0) {
        std::filesystem::path folder_path(folder);
        safe_create_directory(folder_path);
        for (auto& x : std::filesystem::directory_iterator(folder_path)) {
            std::filesystem::remove_all(x);
        }
    }
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

    std::string timing_folder;

    std::string syndrome_folder;
    uint64_t syndrome_rec_hw;

    if (!pp.get_uint32("d", d))             return 1;
    if (!pp.get_float("p", p))              return 1;
    if (!pp.get_uint32("r", r))             return 1;
    if (!pp.get_uint32("blk", blk))         return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_uint64("shots", shots))     return 1;

    bool record_timing_data = pp.get_string("timing-outdir", timing_folder);
    bool record_syndrome_data = pp.get_string("syndrome-outdir", syndrome_folder)
                                    && pp.get_uint64("syndrome-record-hw", syndrome_rec_hw);

    // Define circuits.
    stim::CircuitGenParameters circ_params(r, d, "rotated_memory_z");
    set_error_rate_to(p, circ_params);
    stim::Circuit circuit = stim::generate_surface_code_circuit(circ_params).circuit;

    AstreaDecoder<PyMatching, 8> astrea(circuit);
    PyMatching pymatching(circuit);

    Decoder* base_dec = &astrea;
    if (pp.option_set("pym"))   base_dec = &pymatching;

    BlockDecoder dec(circuit, base_dec, blk);
    if (pp.option_set("pym"))   dec.config.auto_to_global_decoder = true;

    // Add output streams to collect timing data if necessary.
    if (record_timing_data) {
        // Remove folder if it exists -- we want a fresh folder.
        build_folder(timing_folder);
        MPI_Barrier(MPI_COMM_WORLD);
        std::string filename = "proc_" + std::to_string(world_rank) + ".bin";

        dec.config.timing_io.record_data = true;
        dec.config.timing_io.fout = std::ofstream(timing_folder + "/" + filename);
    }

    if (record_syndrome_data) {
        build_folder(syndrome_folder);
        MPI_Barrier(MPI_COMM_WORLD);
        std::string filename = "proc_" + std::to_string(world_rank) + ".txt";

        dec.config.syndrome_io.record_data = true;
        dec.config.syndrome_io.fout = std::ofstream(syndrome_folder + "/" + filename);

        MWPMDecoder* mwpm_base = new MWPMDecoder(circuit);
        dec.config.syndrome_io.base = mwpm_base;
        dec.config.syndrome_io.hw_trigger = syndrome_rec_hw;
    }

    // Setup experiment.
    experiments::G_SHOTS_PER_BATCH = 1'000'000;
    experiments::G_FILTER_OUT_SYNDROMES = false;
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
    fp_t total_ratio_mbhw_thw;
    fp_t total_ratio_mbhw_thw_sqr;
    fp_t total_ratio_mbhw_thw_max;
    uint64_t total_shots_above_blk_th;

    uint64_t total_hw_in_block;
    uint64_t total_hw_sqr_in_block;
    uint64_t max_hw_in_block;
    uint64_t total_blocks;

    uint64_t total_local_uses;
    uint64_t total_global_uses;
    uint64_t total_local_uses_sqr;
    uint64_t total_global_uses_sqr;
    uint64_t max_local_uses;
    uint64_t max_global_uses;

    MPI_Reduce(&dec.total_ratio_mbhw_thw, &total_ratio_mbhw_thw, 1,
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_ratio_mbhw_thw_sqr, &total_ratio_mbhw_thw_sqr, 1,
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_ratio_mbhw_thw_max, &total_ratio_mbhw_thw_max, 1,
            MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_shots_above_blk_th, &total_shots_above_blk_th, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&dec.total_hw_in_block, &total_hw_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_hw_sqr_in_block, &total_hw_sqr_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.max_hw_in_block, &max_hw_in_block, 1,
            MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_blocks, &total_blocks, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&dec.total_local_uses, &total_local_uses, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_global_uses, &total_global_uses, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_local_uses_sqr, &total_local_uses_sqr, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.total_global_uses_sqr, &total_global_uses_sqr, 1,
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.max_local_uses, &max_local_uses, 1,
            MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dec.max_global_uses, &max_global_uses, 1,
            MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        fp_t mean_ratio_mbhw_thw = BMEM_MEAN(total_ratio_mbhw_thw, total_shots_above_blk_th);
        fp_t std_ratio_mbhw_thw = BMEM_STD(total_ratio_mbhw_thw,
                                        total_ratio_mbhw_thw_sqr,
                                        total_shots_above_blk_th);

        fp_t mean_blkhw = BMEM_MEAN(total_hw_in_block, total_blocks);
        fp_t std_blkhw = BMEM_STD(total_hw_in_block,
                                    total_hw_sqr_in_block,
                                    total_blocks);
        fp_t mean_local_uses = BMEM_MEAN(total_local_uses, shots);
        fp_t mean_global_uses = BMEM_MEAN(total_global_uses, shots);
        fp_t std_local_uses = BMEM_STD(total_local_uses, total_local_uses_sqr, shots);
        fp_t std_global_uses = BMEM_STD(total_global_uses, total_global_uses_sqr, shots);

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
                    << "Max Block HW to Total HW Ratio Mean,"
                    << "Max Block HW to Total HW Ratio Std,"
                    << "Max Block HW to Total HW Ratio Max,"
                    << "Block Hamming Weight Mean,"
                    << "Block Hamming Weight Std,"
                    << "Block Hamming Weight Max,"
                    << "Mean Local Decoder Accesses Per Syndrome,"
                    << "Std Local Decoder Accesses Per Syndrome,"
                    << "Max Local Decoder Accesses Per Syndrome,"
                    << "Mean Global Decoder Accesses Per Syndrome,"
                    << "Std Global Decoder Accesses Per Syndrome,"
                    << "Max Global Decoder Accesses Per Syndrome\n";
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
            << mean_ratio_mbhw_thw << ","
            << std_ratio_mbhw_thw << ","
            << total_ratio_mbhw_thw_max << ","
            << mean_blkhw << ","
            << std_blkhw << ","
            << max_hw_in_block << ","
            << mean_local_uses << ","
            << std_local_uses << ","
            << max_local_uses << ","
            << mean_global_uses << ","
            << std_global_uses << ","
            << max_global_uses << "\n";
    }
    
    if (dec.config.syndrome_io.record_data) {
        delete dec.config.syndrome_io.base;
    }

    MPI_Finalize();
}
