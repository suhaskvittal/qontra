/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 *
 *  Memory experiment main file.
 * */

#include "defs.h"
#include "parsing/cmd.h"

#include <stim.h>

#include <mpi.h>

#define SAFE_WR(buf)        if (world_rank == 0) { buf
#define SAFE_WR_FIN         }
#define MPI_DECL(typ, x)    typ x; typ __x
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

    std::string help = 
    "Arguments list:\n"
    + "\tMemory experiment type ( X (-x) Z (-z) default is memory z)\n"
    + "\tDistance (--distance, 32-bit)\n"
    + "\tRounds (--rounds, 32-bit, default is code distance)\n"
    + "\tPhysical error rate mean (--error-rate, float) -- lognormal distribution\n"
    + "\tPhysical error rate standard deviation (--error-rate-std, float)\n"
    + "\tStim file (--stim-file, string) (overrides --distance and --error-rate)\n"
    + "\tShots (--shots, 64-bit)\n"
    + "\tShots per batch (--shots-per-batch, 64-bit)\n"
    + "\tDecoder (MWPM, ...) (--decoder, string)\n"
    + "\tOutput file (--output-file)\n"
    + "\tOutput file type (-csv or -tsv, default csv)\n"
    + "\tSeed (--seed, default 0)\n";

    bool help_requested = option_set("h");
    if (help_requested) {
help_exit:
        SAFE_WR(std::cout) << help; SAFE_WR_FIN
        return 0;
    }

    uint32_t distance;
    uint32_t rounds;
    fp_t pmean;
    fp_t pstd;
    std::string stim_file;
    uint64_t shots;
    uint64_t shots_per_batch;
    std::string decoder;
    std::string output_file;
    uint64_t seed;

    bool is_memory_x = option_set("x");
    bool is_tsv =      option_set("tsv");

    stim::Circuit circuit;
    if (get_string("stim-file", stim_file)) {
        FILE* fptr = fopen(stim_file.c_str());
        circuit = stim::Circuit.from_file(fptr);
        // Null out the other params
        distance = 0;
        pmean = 0;
        pstd = 0;
        fclose(fptr);
    } else if (get_uint32("distance", distance) && get_float("error-rate", pmean)) {
        if (!get_uint32("rounds", rounds)) {
            rounds = distance;
        }
        std::string task = is_memory_x ? "rotated_memory_x" : "rotated_memory_z";
        stim::CircuitGenParams params(rounds, distance, task);
        if (get_float("error-rate-std", pstd)) {
            params.after_clifford_depolarization_stddev = pstd;
            params.before_round_data_depolarization_stddev = pstd;
            params.before_measure_flip_probability_stddev = pstd;
            params.after_reset_flip_probability_stddev = pstd;
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

    if (!get_uint64("shots", shots)) goto help_exit;
    if (!get_uint64("shots-per-batch", shots_per_batch)) {
        shots_per_batch = 100'000;
    }
    if (!get_string("decoder", decoder)) {
        decoder = "MWPM";
    }
    if (!get_string("output-file", output_file)) goto help_exit;
    if (!get_uint64("seed", seed)) seed = 0;

    std::filesystem::path output_path(output_file);
    std::filesystem::path output_folder = output_file.parent_path();
    safe_create_directory(output_folder);
    
    bool write_header = !std::filesystem::exists(output_path);
    std::ofstream out(output_path, std::ios::app);
    
    std::string br(is_tsv ? "\t" : ",");

    if (write_header) { // Create header for stats
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
            << "Max Hamming Weight" << "\n"
            << "Mean Time" << br
            << "Time Std" << br
            << "Max Time" << "\n";
    }
    
    SAFE_WR(out) << distance << br
        << rounds << br
        << pmean << br
        << pstd << br
        << stim_file << br
        << shots << br
        << decoder; SAFE_WR_FIN
    
    // Build the decoder
    decoder::Decoder* dec;
    if (decoder == "MWPM" || decoder == "mwpm") {
        dec = new decoder::MWPMDecoder(circuit);
    }

    // Start the experiment.
    //
    // This code uses MPI to distribute the task.
    uint64_t __shots = shots / world_size;
    if (world_rank == 0)    __shots += shots % world_size;

    std::mt19937_64 rng(seed + world_rank);
    // Local statistics
    MPI_DECL(uint64_t, logical_errors);
    MPI_DECL(uint64_t, hw_sum);
    MPI_DECL(uint64_t, hw_sqr_sum);
    MPI_DECL(uint64_t, hw_max);
    MPI_DECL(fp_t, t_sum);
    MPI_DECL(fp_t, t_sqr_sum);
    MPI_DECL(fp_t, t_max);

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
            stim::simd_bits_range_ref subrow = row.prefix_ref(n_detectors);
            const uint hw = subrow.popcnt();
            if (hw <= 2) continue;
            auto syndrome = decoder::syndrome_to_vector(row, sample_width);
            auto res = dec->decode_error(syndrome);
            // Update stats
            __logical_errors += res.is_error;
            __hw_sum += hw.popcnt();
            __hw_sqr_sum += SQR(hw);
            __hw_max = hw > __hw_max ? hw : __hw_max;
            __t_sum += res.exec_time;
            __t_sqr_sum += SQR(res.exec_time);
            __t_max = res.exec_time > __t_max : res.exec_time : __t_max;
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
    MPI_Reduce(&__t_sum, &t_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__t_sqr_sum, &t_sqr_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__t_max, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // Calculate means and std
    fp_t hw_mean = MEAN(hw_sum, shots);
    fp_t hw_std = STD(hw_mean, hw_sqr_sum, shots);
    fp_t t_mean = MEAN(t_sum, shots);
    fp_t t_std = STD(t_mean, t_sqr_sum, shots);
    // Write the data
    SAFE_WR(out) << hw_mean << br
        << hw_std << br
        << hw_max << br
        << t_mean << br
        << t_std << br
        << t_max << "\n";   SAFE_WR_FIN
    delete dec;

    MPI_Finalize();
}
