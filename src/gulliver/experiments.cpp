/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "gulliver/experiments.h"

std::filesystem::path data_folder(std::string(HOME_DIRECTORY) + "/data");

const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 GULLIVER_RNG(seed);

static const fp_t DEFAULT_ERROR_MEAN = 1e-4;
static const fp_t DEFAULT_ERROR_STDDEV = 30e-5;
static const uint32_t DEFAULT_SHOTS = 100000;

static const uint MQSIM_NQUBITS = 10;
static const uint MQSIM_NDECODER = 5;
static const uint MQSIM_DISTANCE = 5;

static GulliverParams GULLIVER_DEFAULT = {
    4,      // bfu_fetch_width
    10,     // bfu_hw_threshold
    // Memory parameters
    128,    // Number of registers
    // DRAM parameters
    std::string(HOME_DIRECTORY) + "/dramsim3/configs/DDR4_4Gb_x16_1866.ini",
    std::string(HOME_DIRECTORY) + "/src/gulliver/logs",
    1e9,      // clock frequency
};

void
decoder_sram_experiment() {
    std::cout << "Running SRAM experiment...\n";
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    // We examine this for MWPM (motivation).
    for (uint code_dist = 3; code_dist <= 21; code_dist += 2) {
        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, r);
        // MWPM
        MWPMDecoder mwpm_decoder(surf_code_circ);
        // Gulliver
        GulliverParams gulliver_params = GULLIVER_DEFAULT;
        Gulliver gulliver_decoder(surf_code_circ, code_dist*code_dist-1,
                gulliver_params);

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
        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, r);
        // Build decoders.
        // MWPM
        MWPMDecoder mwpm_decoder(surf_code_circ);
        // Gulliver
        GulliverParams gulliver_params = GULLIVER_DEFAULT;
        Gulliver gulliver_decoder(surf_code_circ, code_dist*code_dist-1,
                gulliver_params);
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
        std::cout << "\t\t" << gulliver_decoder.name() << " accessed MWPM " 
            << gulliver_decoder.n_mwpm_accesses << " times out of "
            << gulliver_decoder.n_total_accesses << ".\n";
        std::cout << "\t\t" << "Max cycles in prefetch: " 
            << gulliver_decoder.max_prefetch_cycles << ".\n";
        std::cout << "\t\t" << "Max cycles in BFU: "
            << gulliver_decoder.max_bfu_cycles << ".\n";
        std::cout << "\t\t" << "Max latency: "
            << gulliver_decoder.max_latency << "ns.\n";
        std::cout << "\t\t" << "Number of Rowhammer flips: "
            << gulliver_decoder.simulator->rowhammer_flips() << ".\n";
        std::cout << "\tAdditional stats:\n";
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
    fp_t r = DEFAULT_ERROR_STDDEV;
    for (uint code_dist = 3; code_dist <= 13; code_dist += 2) {
        std::filesystem::path filename(
                "gulliver_timing_d=" + std::to_string(code_dist) + ".txt");
        std::filesystem::path timing_analysis_file = data_folder/filename;

        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, r);
        // Setup decoder.
        GulliverParams decoder_params = GULLIVER_DEFAULT;

        Gulliver decoder(surf_code_circ, code_dist*code_dist-1, decoder_params);
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
    fp_t r = DEFAULT_ERROR_STDDEV;
    for (uint code_dist = 3; code_dist <= 9; code_dist += 2) {
        std::filesystem::path filename(
                "mwpm_timing_d=" + std::to_string(code_dist) + ".txt");
        std::filesystem::path timing_analysis_file = data_folder/filename;

        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, r);
        // Setup decoder.
        MWPMDecoder decoder(surf_code_circ);
        _timing_analysis(timing_analysis_file, &decoder, DEFAULT_SHOTS);
    }
}

void
gulliver_multiple_qubits_experiment() {
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;

    std::vector<stim::Circuit> circuits;
    for (uint i = 0; i < MQSIM_NQUBITS; i++) {
        circuits.push_back(_make_surface_code_circuit(MQSIM_DISTANCE, p, r));
    }

    GulliverParams params = GULLIVER_DEFAULT;
    GulliverMultiQubitSimulator * mqsim = 
        new GulliverMultiQubitSimulator(circuits, MQSIM_NDECODER,
                MQSIM_DISTANCE*MQSIM_DISTANCE-1, GULLIVER_DEFAULT);
    mqsim->benchmark(100000, GULLIVER_RNG);

    std::cout << "Gulliver on " << MQSIM_NQUBITS << " logical qubits across "
        << MQSIM_NDECODER << " decoders for d=" << MQSIM_DISTANCE << " had " 
        << mqsim->n_timeouts << " timeouts, " 
        << mqsim->n_overflows << " overflows, "
        << mqsim->n_uncomputable 
        << " syndromes that could not be computed on the bfu,"
        << " and a max latency of "
        << mqsim->max_latency << "ns.\n";
    delete mqsim;
}

void
mwpm_nonuniform_error_experiment() {
    uint32_t shots = 1000000000;  // 1 billion 
    uint32_t shots_per_round = DEFAULT_SHOTS;
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = DEFAULT_ERROR_STDDEV;
    for (uint code_dist = 3; code_dist <= 11; code_dist += 2) {
        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, r);
        auto proxy_circ = _make_surface_code_circuit(code_dist, p, 0);

        uint32_t nle_proxy = 0;
        uint32_t nle_actual = 0;
#ifdef USE_OMP
#pragma omp parallel for reduction (+:nle_proxy, nle_actual)
#endif
        for (uint32_t batch = 0; batch < shots/shots_per_round; batch++) {
            MWPMDecoder decoder_proxy(proxy_circ);
            MWPMDecoder decoder_actual(proxy_circ);
            // Decoder should've initialized everything for the proxy_circ.
            // We replace the circuit with the actual circuit.
            decoder_actual.circuit = surf_code_circ;
            b_decoder_ler(&decoder_proxy, shots_per_round, GULLIVER_RNG, false);
            b_decoder_ler(&decoder_actual, shots_per_round, GULLIVER_RNG, false);
            nle_proxy += decoder_proxy.n_logical_errors;
            nle_actual += decoder_actual.n_logical_errors;
        }
        fp_t ler_proxy = ((fp_t) nle_proxy) / ((fp_t) shots);
        fp_t ler_actual = ((fp_t) nle_actual) / ((fp_t) shots);
        std::cout << "[d = " << code_dist << "] MWPM LER = " << ler_actual
            << " (" << ler_proxy << ").\n";
    }
}

void
surface_code_hamming_weight_experiment() {
//  uint32_t shots = 1000000000; // one billion.
    uint32_t shots = 100000000;
    uint32_t shots_per_round = DEFAULT_SHOTS;
    fp_t p = DEFAULT_ERROR_MEAN;
    fp_t r = 0;DEFAULT_ERROR_STDDEV;

    for (uint code_dist = 3; code_dist <= 11; code_dist += 2) {
        std::cout << "Distance " << code_dist << ":\n";

        std::map<uint, uint32_t> hwfreq_table;
        std::map<uint, fp_t> maxprob_table;

        auto surf_code_circ = _make_surface_code_circuit(code_dist, p, 0);
        uint n_detectors = surf_code_circ.count_detectors();
        uint n_observables = surf_code_circ.count_observables();
        MWPMDecoder decoder(surf_code_circ);
    
        b_decoder_ler(&decoder, shots, GULLIVER_RNG, false);
        fp_t ler = ((fp_t)decoder.n_logical_errors) / ((fp_t)shots);
        std::cout << "\tMWPM LER = " << ler << "\n";
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (uint32_t batch = 0; batch < shots/shots_per_round; batch++) {
            stim::simd_bit_table sample_buffer
                = stim::detector_samples(surf_code_circ, shots_per_round, 
                        false, true, GULLIVER_RNG);
            sample_buffer = sample_buffer.transposed();
            for (uint32_t sn = 0; sn < shots_per_round; sn++) {
                std::vector<uint8_t> syndrome = 
                    _to_vector(sample_buffer[sn], n_detectors, n_observables);
                uint hw = std::accumulate(
                        syndrome.begin(),
                        syndrome.begin() + n_detectors,
                        0);
                DecoderShotResult res = decoder.decode_error(syndrome);
#ifdef USE_OMP
#pragma omp critical
#endif
                {
                    if (hwfreq_table.count(hw) == 0) {
                        hwfreq_table[hw] = 0;
                        maxprob_table[hw] = std::numeric_limits<fp_t>::max();
                    }
                    hwfreq_table[hw]++;

                    fp_t matching_weight = 0.0;
                    for (auto assign : res.matching) {
                        matching_weight += decoder.path_table[assign].distance;
                    }
                    if (matching_weight < maxprob_table[hw]) {
                        maxprob_table[hw] = matching_weight*0.5;
                    }
                }
            }
        }
        for (auto kv_entry : hwfreq_table) {
            uint hw = kv_entry.first;
            uint32_t f = kv_entry.second;
            std::cout << "\tHamming weight " << hw << " has frequency "
                << f << " (max probability = " 
                << 1.0/(pow(M_E, maxprob_table[hw])+1) << ")\n";
        }
    }
}

stim::Circuit
_make_surface_code_circuit(uint code_dist, fp_t mean, fp_t stddev) {
    stim::CircuitGenParameters surf_code_params(
            code_dist, code_dist, "rotated_memory_z");
    // Setup error rates.
    surf_code_params.after_clifford_depolarization = mean;
    surf_code_params.before_round_data_depolarization = mean;
    surf_code_params.before_measure_flip_probability = mean;
    surf_code_params.after_reset_flip_probability = mean;

    surf_code_params.after_clifford_depolarization_stddev = stddev;
    surf_code_params.before_round_data_depolarization_stddev = stddev;
    surf_code_params.before_measure_flip_probability_stddev = stddev;
    surf_code_params.after_reset_flip_probability_stddev = stddev;

    stim::Circuit circ = generate_surface_code_circuit(surf_code_params).circuit;

    return circ;
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

