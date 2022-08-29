/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "gulliver/experiments.h"

std::filesystem::path data_folder(std::string(HOME_DIRECTORY) + "/data");

const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 GULLIVER_RNG(seed);

static fp_t DEFAULT_ERROR_MEAN = 1e-4;
static fp_t DEFAULT_ERROR_STDDEV = 15e-5;
static uint32_t DEFAULT_SHOTS = 100000;

static GulliverParams GULLIVER_DEFAULT = {
    1,      // n_bfu
    5,      // n_bfu_cycles_per_add
    7,      // bfu_hw_threshold
    500e6,  // clock_frequency
    // Cache parameters,
    14, // C, cache size is 2**C
    2,  // S
    10,  // B
    9,  // tlb C
    4,  // tlb B
    // DRAM parameters
    std::string(HOME_DIRECTORY) + "/dramsim3/configs/DDR3_1Gb_x8_1333.ini",
    std::string(HOME_DIRECTORY) + "/src/gulliver/logs"
};

void
decoder_sram_experiment() {
    std::cout << "Running SRAM experiment...\n";
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    // We examine this for MWPM (motivation).
    for (uint code_dist = 3; code_dist <= 21; code_dist += 2) {
        stim::CircuitGenParameters surf_code_params(
                code_dist, code_dist, "rotated_memory_z");
        // Setup error rates.
        surf_code_params.after_clifford_depolarization = p;
        surf_code_params.before_round_data_depolarization = p;
        surf_code_params.before_measure_flip_probability = p;
        surf_code_params.after_reset_flip_probability = p;

        surf_code_params.after_clifford_depolarization_stddev = r;
        surf_code_params.before_round_data_depolarization_stddev = r;
        surf_code_params.before_measure_flip_probability_stddev = r;
        surf_code_params.after_reset_flip_probability_stddev = r;

        stim::Circuit surf_code_circ = 
            generate_surface_code_circuit(surf_code_params).circuit;
        // MWPM
        MWPMDecoder mwpm_decoder(surf_code_circ);
        // Gulliver
        GulliverParams gulliver_params = GULLIVER_DEFAULT;
        Gulliver gulliver_decoder(surf_code_circ, gulliver_params);

        Decoder * decoder_array[] = {
            &mwpm_decoder, 
            &gulliver_decoder
        };

        for (Decoder * d_p : decoder_array) {
            uint64_t n_bytes_sram = d_p->sram_cost();
            std::cout << "For distance " << code_dist 
                << ", " << d_p->name() << " requires "
                << n_bytes_sram << " bytes of SRAM, or " 
                << n_bytes_sram / 1024.0 << " KB, or " 
                << n_bytes_sram / (1024.0*1024.0) << " MB.\n";
        }
    }
}

void
decoder_analysis_experiment() {
    std::cout << "Running decoder analysis experiment...\n";
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    uint32_t shots = 1000000;
    // Setup circuit.
    for (uint code_dist = 3; code_dist <= 11; code_dist += 2) {
        stim::CircuitGenParameters surf_code_params(
                code_dist, code_dist, "rotated_memory_z");
        // Setup error rates.
        surf_code_params.after_clifford_depolarization = p;
        surf_code_params.before_round_data_depolarization = p;
        surf_code_params.before_measure_flip_probability = p;
        surf_code_params.after_reset_flip_probability = p;

        surf_code_params.after_clifford_depolarization_stddev = r;
        surf_code_params.before_round_data_depolarization_stddev = r;
        surf_code_params.before_measure_flip_probability_stddev = r;
        surf_code_params.after_reset_flip_probability_stddev = r;

        stim::Circuit surf_code_circ = 
            generate_surface_code_circuit(surf_code_params).circuit;
        // Build decoders.
        // MWPM
        MWPMDecoder mwpm_decoder(surf_code_circ);
        // Gulliver
        GulliverParams gulliver_params = GULLIVER_DEFAULT;
        Gulliver gulliver_decoder(surf_code_circ, gulliver_params);
        // CliqueDecoder
        CliqueParams clique_params = {
            1,
            1,
            1,
            5.0e9,
            code_dist*code_dist - 1
        };
        CliqueDecoder clique_decoder(surf_code_circ, clique_params); 

        Decoder * decoder_array[] = {
            &mwpm_decoder, 
            &gulliver_decoder,
            &clique_decoder
        };
        std::cout << "Physical error rate: " << p 
            << ", code distance: " << code_dist 
            << ", syndrome size: " << surf_code_circ.count_detectors() 
            <<"\n";
        for (Decoder * d_p : decoder_array) {
            b_decoder_ler(d_p, shots, GULLIVER_RNG, false); 
            fp_t ler = ((fp_t)d_p->n_logical_errors) / ((fp_t)shots);
            std::cout << "\t" << d_p->name() << " LER = " << ler 
                << ", time taken: " << d_p->mean_execution_time 
                << " (max time = " << d_p->max_execution_time 
                << "; " << d_p->max_execution_time_for_correctable
                << ")\n";
        }
        std::cout << "\tGulliver stats:\n";
        fp_t cache_miss_rate = 
            ((fp_t) gulliver_decoder.cache->n_misses) 
            / gulliver_decoder.cache->n_accesses;
        std::cout << "\t\tCache misses: " << gulliver_decoder.cache->n_misses
            << " from " << gulliver_decoder.cache->n_accesses
            << " accesses (miss rate = " << cache_miss_rate 
            << ").\n";
        fp_t tlb_miss_rate = 
            ((fp_t) gulliver_decoder.cache->n_tlb_misses) 
            / gulliver_decoder.cache->n_tlb_accesses;
        std::cout << "\t\tTLB misses: " << gulliver_decoder.cache->n_tlb_misses
            << " from " << gulliver_decoder.cache->n_tlb_accesses
            << " accesses (miss rate = " << tlb_miss_rate 
            << ").\n";
        std::cout << "\tAdditional stats:\n";
        std::cout << "\t\t" << gulliver_decoder.name() << " accessed MWPM " 
            << gulliver_decoder.n_mwpm_accesses << " times out of "
            << gulliver_decoder.n_total_accesses << ".\n";
        std::cout << "\t\t" << clique_decoder.name() << " accessed MWPM "
            << clique_decoder.n_mwpm_accesses << " times out of "
            << clique_decoder.n_total_accesses << ".\n";
    }
}

void
gulliver_timing_experiment() {
    std::cout << "Running Gulliver timing analysis experiment...\n";

    create_directory(data_folder);

    fp_t p = DEFAULT_ERROR_MEAN;
    for (uint code_dist = 3; code_dist <= 13; code_dist += 2) {
        std::filesystem::path filename(
                "gulliver_timing_d=" + std::to_string(code_dist) + ".txt");
        std::filesystem::path timing_analysis_file = data_folder/filename;
        // Setup Stim parameters.
        stim::CircuitGenParameters surf_code_params(
                code_dist, code_dist, "rotated_memory_z");
        // Setup error rates.
        surf_code_params.after_clifford_depolarization = p;
        surf_code_params.before_round_data_depolarization = p;
        surf_code_params.before_measure_flip_probability = p;
        surf_code_params.after_reset_flip_probability = p;

        stim::Circuit surf_code_circ = 
            generate_surface_code_circuit(surf_code_params).circuit;
        // Setup decoder.
        GulliverParams decoder_params = GULLIVER_DEFAULT;

        Gulliver decoder(surf_code_circ, decoder_params);
        _timing_analysis(timing_analysis_file, &decoder, DEFAULT_SHOTS);

        uint32_t n_nonzero_syndromes = 0;
        auto completable_syndromes = 
            decoder.nonzero_syndromes_completed_within(1000, n_nonzero_syndromes);
        std::cout << "Distance = " << code_dist
            << "\tGulliver can complete "
            << completable_syndromes.size() << " of "
            << n_nonzero_syndromes << " within 1us.\n";
    }
}

void
mwpm_timing_experiment() {
    std::cout << "Running MWPM timing analysis experiment...\n";

    create_directory(data_folder);

    fp_t p = DEFAULT_ERROR_MEAN;
    for (uint code_dist = 3; code_dist <= 9; code_dist += 2) {
        std::filesystem::path filename(
                "mwpm_timing_d=" + std::to_string(code_dist) + ".txt");
        std::filesystem::path timing_analysis_file = data_folder/filename;
        // Setup Stim parameters.
        stim::CircuitGenParameters surf_code_params(
                code_dist, code_dist, "rotated_memory_z");
        // Setup error rates.
        surf_code_params.after_clifford_depolarization = p;
        surf_code_params.before_round_data_depolarization = p;
        surf_code_params.before_measure_flip_probability = p;
        surf_code_params.after_reset_flip_probability = p;

        stim::Circuit surf_code_circ = 
            generate_surface_code_circuit(surf_code_params).circuit;
        // Setup decoder.
        MWPMDecoder decoder(surf_code_circ);
        _timing_analysis(timing_analysis_file, &decoder, DEFAULT_SHOTS);
    }
}

void
mwpm_sweep_experiment() {
    std::cout << "Running sweep analysis experiment...\n";

    std::filesystem::path folder("mwpm_sweep");
    std::filesystem::path sweep_analysis_folder = data_folder/folder;
    // Create folder if it does not exist.
    create_directory(data_folder);
    create_directory(sweep_analysis_folder);

    ErrorThresholdSweepParams::decoder_gen_f dgf = 
        [] (uint code_dist, fp_t error) {
            stim::CircuitGenParameters circ_params(
                    code_dist, code_dist, "rotated_memory_z");
            // Setup error rates.
            circ_params.after_clifford_depolarization = error;
            circ_params.before_round_data_depolarization = error;
            circ_params.before_measure_flip_probability = error;
            circ_params.after_reset_flip_probability = error;
            stim::Circuit circ = 
                generate_surface_code_circuit(circ_params).circuit;
            Decoder * decoder = new MWPMDecoder(circ);
            return decoder;
        };

    ErrorThresholdSweepParams params = {
        dgf,
        1e-3,
        1e-2,
        1e-3,
        3,
        11,
        2
    };
    _threshold_sweep(sweep_analysis_folder, params, DEFAULT_SHOTS);
}

void
gulliver_cache_experiment() {
    std::filesystem::path folder("gulliver_cache");
    std::filesystem::path cache_analysis_folder(data_folder/folder);
    // Create directories if they do not exist.
    create_directory(data_folder);
    create_directory(cache_analysis_folder);

    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    uint32_t shots = 100000;
    for (uint code_dist = 3; code_dist <= 7; code_dist += 2) {
        std::filesystem::path out_file(
                "dist" + std::to_string(code_dist)  + ".txt");
        std::filesystem::path output_path = cache_analysis_folder/out_file;
        // Build circuit
        stim::CircuitGenParameters surf_code_params(
                code_dist, code_dist, "rotated_memory_z");
        // Setup error rates.
        surf_code_params.after_clifford_depolarization = p;
        surf_code_params.before_round_data_depolarization = p;
        surf_code_params.before_measure_flip_probability = p;
        surf_code_params.after_reset_flip_probability = p;

        surf_code_params.after_clifford_depolarization_stddev = r;
        surf_code_params.before_round_data_depolarization_stddev = r;
        surf_code_params.before_measure_flip_probability_stddev = r;
        surf_code_params.after_reset_flip_probability_stddev = r;

        stim::Circuit surf_code_circ = 
            generate_surface_code_circuit(surf_code_params).circuit;
        // Now, sweep over C, B, S.
        std::cout << "Distance " << code_dist << ":\n";
        _cache_sweep(output_path, surf_code_circ, shots);
    }
}

void 
_timing_analysis(const std::filesystem::path& output_file,
        Decoder * decoder_p, uint32_t shots)
{
    // Benchmark decoder.
    b_decoder_ler(decoder_p, shots, GULLIVER_RNG);
    auto syndromes = decoder_p->nonzero_syndromes_and_time_taken();
    // Write all fast syndromes to file.
    // Format: <syndrome> <time_taken>
    std::ofstream out(output_file);
    for (auto pair : syndromes) {
        std::vector<uint8_t> syndrome = pair.first;
        fp_t time_taken = pair.second;
        std::string bitvec;
        for (uint8_t bit : syndrome) {
            if (bit) {
                bitvec.push_back('1');
            } else {
                bitvec.push_back('0');
            }
        }
        out << bitvec << " " << time_taken << "\n";
    }
}

void 
_threshold_sweep(const std::filesystem::path& output_folder,
        const ErrorThresholdSweepParams& params, uint32_t shots) 
{
    auto data = sweep_error_threshold(params, shots, GULLIVER_RNG);
    // Write data to output folder.
    std::filesystem::path circ_output_folder = output_folder/"circuits";
    std::filesystem::path sweep_output_file = output_folder/"sweep.txt";
    // Create any folders that do not exist.
    safe_create_directory(output_folder);
    safe_create_directory(circ_output_folder);
    // Finally, run the sweep.
    ErrorThresholdOutputParams output_params = {
        true,
        circ_output_folder
    };
    std::ofstream out(sweep_output_file);
    write_sweep_data(data, out, output_params);
}

void
_cache_sweep(const std::filesystem::path& output_file,
        const stim::Circuit& circuit, uint32_t shots)
{
    std::ofstream out(output_file);
    for (uint C = 10; C <= 20; C++) {
        for (uint B = 6; B <= 12; B++) {
            for (uint S = 0; S <= C - B; S++) {
                GulliverParams params = GULLIVER_DEFAULT;
                params.cacheC = C;
                params.cacheB = B;
                params.cacheS = S;
                Gulliver decoder(circuit, params);
                b_decoder_ler(&decoder, shots, GULLIVER_RNG, false); 
                fp_t cache_miss_rate = 
                    ((fp_t) decoder.cache->n_misses) 
                    / decoder.cache->n_accesses;
                out << C << " " << S << " " << B 
                    << " " << cache_miss_rate << "\n";
                std::cout << "\t" << C << " " << S << " " << B 
                    << " " << cache_miss_rate << "\n";
            }
        }
    }
}
