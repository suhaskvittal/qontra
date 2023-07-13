/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 *
 *  Memory experiment main file.
 * */

#include "decoder/decoder.h"
#include "decoder/mwpm.h"
#include "decoder/window.h"
#include "defs.h"
#include "experiments.h"
#include "instruction.h"
#include "parsing/cmd.h"

#include <stim.h>

#include <fstream>
#include <iostream>

#include <mpi.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tMemory experiment type (X (-x) Z (-z) default is -z)\n";
    help += "\tDistance (--distance)\n";
    help += "\tRounds (--rounds, default is code distance)\n";
    help += "\tPhysical error rate mean (--error-rate)\n";
    help += "\tPhysical error rate standard deviation (--error-rate-std)\n";
    help += "\tStim file (--stim-file, string)";
    help += " (overrides distance and error-rate)\n";
    help += "\tShots (--shots)\n";
    help += "\tShots per batch (--shots-per-batch)\n";
    help += "\tDecoder (MWPM, ...) (--decoder)\n";
    help += "\tOutput file (--output-file)\n";
    help += "\tOutput file type (-csv or -tsv, default csv)\n";
    help += "\tUsing sliding window (-window)\n";
    help += "\tSliding window stim file (--wstim-file,";
    help += " must provide if --stim-file was used, overrides --wr)\n";
    help += "\tSliding window rounds per decode (--wr)\n";
    help += "\tSliding window detectors per round (--wdpr)\n";
    help += "\tSliding window commit window (--wcomm)\n";
    help += "\tSeed (--seed, default 0)\n";

    bool help_requested = parser.option_set("h");
    if (help_requested) {
help_exit:
        if (world_rank == 0) {
            std::cout << help;
        }
        return 0;
    }

    uint32_t distance;
    uint32_t rounds = 0;
    fp_t pmean;
    fp_t pstd = 0;
    std::string stim_file;
    uint64_t shots;
    uint64_t shots_per_batch;
    std::string decoder;
    std::string output_file;
    uint64_t seed;

    bool is_memory_x =  parser.option_set("x");
    bool is_tsv =       parser.option_set("tsv");
    bool use_window =   parser.option_set("window");

    stim::Circuit circuit;

    if (parser.get_string("stim-file", stim_file)) {
        FILE* fptr = fopen(stim_file.c_str(), "r");
        circuit = stim::Circuit::from_file(fptr);
        // Null out the other params
        distance = 0;
        pmean = 0;
        pstd = 0;
        fclose(fptr);
    } else if (parser.get_uint32("distance", distance) 
            && parser.get_float("error-rate", pmean)) 
    {
        if (!parser.get_uint32("rounds", rounds)) {
            rounds = distance;
        }
        std::string task = is_memory_x ? "rotated_memory_x" : "rotated_memory_z";
        stim::CircuitGenParameters params(rounds, distance, task);
        if (parser.get_float("error-rate-std", pstd)) {
            params.after_clifford_depolarization_stddev = pstd;
            params.before_round_data_depolarization_stddev = pstd;
            params.before_measure_flip_probability_stddev = pstd;
            params.after_reset_flip_probability_stddev = pstd;
        } else {
            pstd = 0;
        }
        params.after_clifford_depolarization = pmean;
        params.before_round_data_depolarization = pmean;
        params.before_measure_flip_probability = pmean;
        params.after_reset_flip_probability = pmean;

        circuit = stim::generate_surface_code_circuit(params).circuit;
        // Null out stim_file
        stim_file = "";
    } else {
        goto help_exit;
    }

    uint32_t wd;
    uint32_t wr;
    uint32_t wcomm;
    uint32_t wdpr;
    std::string wfile;
    stim::Circuit wcirc;

    if (use_window) {
        if (parser.get_string("wstim-file", wfile)) {
            FILE* fptr = fopen(stim_file.c_str(), "r");
            wcirc = stim::Circuit::from_file(fptr);
            
            wd = 0;
            wr = 0;
        } else if (parser.get_uint32("distance", wd)
            && parser.get_uint32("wr", wr))
        {
            std::string task = is_memory_x ? "rotated_memory_x" : "rotated_memory_z";
            stim::CircuitGenParameters params(wr+1, wd, task);
            if (parser.get_float("error-rate-std", pstd)) {
                params.after_clifford_depolarization_stddev = pstd;
                params.before_round_data_depolarization_stddev = pstd;
                params.before_measure_flip_probability_stddev = pstd;
                params.after_reset_flip_probability_stddev = pstd;
            } else {
                pstd = 0;
            }
            params.after_clifford_depolarization = pmean;
            params.before_round_data_depolarization = pmean;
            params.before_measure_flip_probability = pmean;
            params.after_reset_flip_probability = pmean;

            wcirc = stim::generate_surface_code_circuit(params).circuit;
        } else {
            goto help_exit;
        }

        if (!parser.get_uint32("wdpr", wdpr))   goto help_exit;
        if (!parser.get_uint32("wcomm", wcomm)) goto help_exit;
    }

    if (!parser.get_uint64("shots", shots)) goto help_exit;
    if (!parser.get_uint64("shots-per-batch", shots_per_batch)) {
        shots_per_batch = 100'000;
    }
    if (!parser.get_string("decoder", decoder)) {
        decoder = "MWPM";
    }
    if (!parser.get_string("output-file", output_file)) goto help_exit;
    if (!parser.get_uint64("seed", seed)) seed = 0;

    std::filesystem::path output_path(output_file);
    std::filesystem::path output_folder = output_path.parent_path();
    safe_create_directory(output_folder);
    
    bool write_header = !std::filesystem::exists(output_path);
    MPI_Barrier(MPI_COMM_WORLD);

    std::ofstream out(output_path, std::ios::app);
    
    std::string br(is_tsv ? "\t" : ",");

    if (write_header && world_rank == 0) { // Create header for stats
        out << "Distance" << br
            << "Rounds" << br
            << "Pmean" << br
            << "Pstd" << br
            << "Stim File" << br
            << "Shots" << br
            << "Decoder" << br
            << "Logical Error Probability" << br
            << "Mean Hamming Weight" << br
            << "Hamming Weight Std" << br
            << "Max Hamming Weight" << br
            << "Mean Time" << br
            << "Time Std" << br
            << "Max Time" << "\n";
    }
    
    if (world_rank == 0) {
        out << distance << br
            << rounds << br
            << pmean << br
            << pstd << br
            << output_path.filename() << br
            << shots << br
            << decoder;
        if (use_window) {
            out << "(r"  << wr << "|c" << wcomm << ")";
        }
        out << br;
    }
    
    // Build the decoder
    if (use_window) std::swap(circuit, wcirc);
    decoder::Decoder* dec;
    decoder::Decoder* base_dec;
    if (decoder == "MWPM" || decoder == "mwpm") {
        dec = new decoder::MWPMDecoder(circuit);
    } else if (decoder == "NONE" || decoder == "none") {
        std::cout << help;
        return 1;
    }
    // If we must decode in a sliding window, then build a sliding window decoder.
    if (use_window) {
        base_dec = dec;
        dec = new decoder::WindowDecoder(wcirc, base_dec, wcomm, wdpr);
    }
    // Execute the experiment
    auto res = memory_experiment(dec, shots);
    // Write the data
    if (world_rank == 0) {
        out << res.logical_error_rate << br
            << res.hw_mean << br
            << res.hw_std << br
            << res.hw_max << br
            << res.t_mean << br
            << res.t_std << br
            << res.t_max << "\n";
    }
    if (use_window) delete base_dec;
    delete dec;
    MPI_Finalize();
    return 0;
}
