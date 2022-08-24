/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "gulliver/experiments.h"

std::filesystem::path data_folder(std::string(HOME_DIRECTORY) + "/data");

const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 GULLIVER_RNG(seed);

static fp_t DEFAULT_ERROR_MEAN = 1e-3;
static fp_t DEFAULT_ERROR_STDDEV = 15e-4;
static uint32_t DEFAULT_SHOTS = 100000;

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
        MWPMDecoder mwpm_decoder(surf_code_circ);
        uint64_t n_bytes_sram = mwpm_decoder.sram_cost();
        std::cout << "For distance " << code_dist 
            << ", a MWPM decoder requires " << n_bytes_sram
            << " bytes of SRAM, or " 
            << n_bytes_sram / 1024.0 << " KB, or " 
            << n_bytes_sram / (1024.0*1024.0) << " MB.\n";
    }
}

void
decoder_analysis_experiment() {
    std::cout << "Running decoder analysis experiment...\n";
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    uint32_t shots = DEFAULT_SHOTS;
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
        GulliverParams gulliver_params = {
            1,      // n_bfu
            5,      // n_bfu_cycles_per_add
            6,      // bfu_hw_threshold
            250e6   // clock_frequency
        };
        GulliverDecoder gulliver_decoder(surf_code_circ, gulliver_params);

        Decoder * decoder_array[] = {&mwpm_decoder, &gulliver_decoder};
        std::cout << "Physical error rate: " << p 
            << ", code distance: " << code_dist << "\n";
        for (Decoder * d_p : decoder_array) {
            b_decoder_ler(d_p, shots, GULLIVER_RNG); 
            fp_t ler = ((fp_t)d_p->n_logical_errors) / ((fp_t)shots);
            std::cout << "\t" << d_p->name() << " LER = " << ler << "\n";
        }
    }
}

void
gulliver_timing_experiment() {
    std::cout << "Running Gulliver timing analysis experiment...\n";

    create_directory(data_folder);

    fp_t p = DEFAULT_ERROR_MEAN;
    for (uint code_dist = 3; code_dist <= 9; code_dist += 2) {
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
        GulliverParams decoder_params = {
            1,      // n_bfu
            5,      // n_bfu_cycles_per_add
            6,      // bfu_hw_threshold
            250e6   // clock_frequency
        };

        GulliverDecoder decoder(surf_code_circ, decoder_params);
        GulliverTimingAnalysisParams params = {
            &decoder,
            100000
        };
        _timing_analysis(timing_analysis_file, params);
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
        GulliverTimingAnalysisParams params = {
            &decoder,
            100000
        };
        _timing_analysis(timing_analysis_file, params);
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

    GulliverSweepAnalysisParams params = {
        dgf,
        1e-3,
        1e-2,
        1e-3,
        3,
        11,
        2,
        100000
    };
    _sweep_analysis(sweep_analysis_folder, params);
}

void 
_timing_analysis(const std::filesystem::path& output_file,
        const GulliverTimingAnalysisParams& params)
{
    // Benchmark decoder.
    b_decoder_ler(params.decoder_p, params.shots, GULLIVER_RNG);
    auto syndromes = params.decoder_p->nonzero_syndromes_and_time_taken();
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
_sweep_analysis(const std::filesystem::path& output_folder,
        const GulliverSweepAnalysisParams& params) 
{
    // Declare sweep params and perform sweep.
    ErrorThresholdSweepParams sweep_params = {
        params.dgf,
        params.error_lb,
        params.error_ub,
        params.error_step,
        params.code_distance_lb,
        params.code_distance_ub,
        params.code_distance_step
    };
    auto data = sweep_error_threshold(sweep_params, params.shots, GULLIVER_RNG);
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
