/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 *
 *  Memory experiment main file.
 * */

#include "decoder/decoder.h"
#include "decoder/mwpm.h"
#include "defs.h"
#include "parsing/cmd.h"

#include <stim.h>

#include <fstream>
#include <iostream>

#include <mpi.h>

#define SQR(x)              (x)*(x)
#define MEAN(s, n)          ((fp_t)(s))/((fp_t)(n))
#define STD(m, ss, n)       ( ((fp_t)(ss))/((fp_t)(n)) - SQR(m) )

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tMemory experiment type (X (-x) Z (-z) default is -z)\n";
    help += "\tDistance (--distance, 32-bit)\n";
    help += "\tRounds (--rounds, 32-bit, default is code distance)\n";
    help += "\tPhysical error rate mean (--error-rate, float, lognormal)\n";
    help += "\tPhysical error rate standard deviation (--error-rate-std, float)\n";
    help += "\tStim file (--stim-file, string) (overrides distance and error-rate)\n";
    help += "\tShots (--shots, 64-bit)\n";
    help += "\tShots per batch (--shots-per-batch, 64-bit)\n";
    help += "\tDecoder (MWPM, ...) (--decoder, string)\n";
    help += "\tOutput file (--output-file)\n";
    help += "\tOutput file type (-csv or -tsv, default csv)\n";
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

    bool is_memory_x = parser.option_set("x");
    bool is_tsv =      parser.option_set("tsv");

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
            << decoder << br;
    }
    
    // Build the decoder
    decoder::Decoder* dec;
    if (decoder == "MWPM" || decoder == "mwpm") {
        dec = new decoder::MWPMDecoder(circuit);
    } else if (decoder == "NONE" || decoder == "none") {
        dec = nullptr;
    }

    // Start the experiment.
    //
    // This code uses MPI to distribute the task.
    uint64_t __shots = shots / world_size;
    if (world_rank == 0)    __shots += shots % world_size;

    std::mt19937_64 rng(seed + world_rank);

    uint64_t    logical_errors, __logical_errors = 0;
    uint64_t    hw_sum, __hw_sum = 0;
    uint64_t    hw_sqr_sum, __hw_sqr_sum = 0;
    uint64_t    hw_max, __hw_max = 0;
    fp_t        t_sum, __t_sum = 0;
    fp_t        t_sqr_sum, __t_sqr_sum = 0;
    fp_t        t_max, __t_max = 0;

    const uint sample_width = 
        circuit.count_detectors() + circuit.count_observables();
    while (__shots > 0) {
        const uint64_t shots_this_batch = 
            __shots < shots_per_batch ? __shots : shots_per_batch;
        auto samples = stim::detector_samples(
                            circuit, shots_this_batch, false, true, rng);
        samples = samples.transposed();
        for (uint64_t t = 0; t < shots_this_batch; t++) {
            stim::simd_bits_range_ref row = samples[t];
            stim::simd_bits_range_ref subrow = 
                row.prefix_ref(circuit.count_detectors());
            const uint hw = subrow.popcnt();

            __hw_sum += hw;
            __hw_sqr_sum += SQR(hw);
            __hw_max = hw > __hw_max ? hw : __hw_max;

            if (hw <= 2 && dec != nullptr) continue;
            auto syndrome = decoder::syndrome_to_vector(row, sample_width);
            auto res = dec->decode_error(syndrome);

            __logical_errors += res.is_error;
            __t_sum += res.exec_time;
            __t_sqr_sum += SQR(res.exec_time);
            __t_max = res.exec_time > __t_max ? res.exec_time : __t_max;
        }
        __shots -= shots_this_batch;
    }

    MPI_Reduce(&__logical_errors, &logical_errors, 1, MPI_UNSIGNED_LONG,
                MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__hw_sum, &hw_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                0, MPI_COMM_WORLD);
    MPI_Reduce(&__hw_sqr_sum, &hw_sqr_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                0, MPI_COMM_WORLD);
    MPI_Reduce(&__hw_max, &hw_max, 1, MPI_UNSIGNED_LONG, MPI_MAX,
                0, MPI_COMM_WORLD);
    MPI_Reduce(&__t_sum, &t_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__t_sqr_sum, &t_sqr_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__t_max, &t_max, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // Calculate means and std
    fp_t ler = ((fp_t)logical_errors) / ((fp_t)shots);
    fp_t hw_mean = MEAN(hw_sum, shots);
    fp_t hw_std = STD(hw_mean, hw_sqr_sum, shots);
    fp_t t_mean = MEAN(t_sum, shots);
    fp_t t_std = STD(t_mean, t_sqr_sum, shots);
    // Write the data
    if (world_rank == 0) {
        std::cout << "errors = " << logical_errors << "\n";
        out << ler << br
            << hw_mean << br
            << hw_std << br
            << hw_max << br
            << t_mean << br
            << t_std << br
            << t_max << "\n";
    }
    if (dec != nullptr) delete dec;

    MPI_Finalize();
}
