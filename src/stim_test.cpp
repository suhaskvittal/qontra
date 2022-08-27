/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "stim_test.h"

const auto seed = std::chrono::system_clock::now()
                            .time_since_epoch().count();
std::mt19937_64 RNG(seed);

void test_circuit_build() {
    std::cout << "======= CIRCUIT BUILD TEST =======\n";
    stim::Circuit circ;
    circ.append_from_text("H 0");
    circ.append_from_text("CNOT 0 1");
    circ.append_from_text("M 0 1");
    // Print out stim circuit.
    std::cout << "Stim circuit:\n";
    std::cout << circ.str() << "\n";
    // Sample circuit.
    stim::simd_bit_table samples = stim::FrameSimulator::sample_flipped_measurements(circ, 10, RNG);
    std::cout << "Sampler results:\n";
    uint8_t prev_results[2] = {0,0};
    std::string det_seq[2];
    for (uint i = 0; i < 10; i++) {
        for (uint64_t j = 0; j <= 1; j++) {
            uint8_t res = samples[j][i] ^ prev_results[j];
            det_seq[j].push_back(res + '0');
            prev_results[j] = res;
        }
    }
    std::cout << "meas sequence 0: " << det_seq[0] << "\n";
    std::cout << "meas sequence 1: " << det_seq[1] << "\n";
}

void test_detector_build() {
    std::cout << "======= DETECTOR BUILD TEST =======\n";
    stim::Circuit circ( // Provide string description.
        "H 0\n"
        "CX 0 1\n"
        "X_ERROR(0.2) 0 1\n"
        "M 0 1\n"
        "DETECTOR rec[-1] rec[-2]"
    );
    std::cout << "Stim circuit:\n";
    std::cout << circ.str() << "\n";
    // Sample detectors.
    stim::simd_bit_table samples = stim::detector_samples(circ, 10, false, false, RNG);
    std::string det_seq;
    for (uint i = 0; i < 10; i++) {
        uint8_t res = samples[0][i];
        det_seq.push_back(res + '0');
    }
    std::cout << "det sequence: " << det_seq << "\n";
}

void test_qec_build() {
    std::cout << "======= QEC BUILD TEST =======\n";
    std::cout << "Testing repetition code...\n";

    stim::CircuitGenParameters rep_code_params(100, 9, "memory");
    rep_code_params.before_round_data_depolarization = 0.03;

    stim::Circuit rep_code_circ = 
        generate_rep_code_circuit(rep_code_params).circuit;
    stim::simd_bit_table rep_code_samples = 
        stim::detector_samples(rep_code_circ, 1, false, false, RNG);
    uint8_t rep_code_n_det = 8;  // As code distance is 9
    for (uint i = 0; i < 100; i++) { 
        std::string detection_event;
        for (uint8_t j = 0; j < rep_code_n_det; j++) {
            if (rep_code_samples[i*rep_code_n_det + j][0]) {
                detection_event.push_back('!');
            } else {
                detection_event.push_back('_');
            }
        }
        std::cout << "Detection event: " << detection_event << "\n";
    }

    std::cout << "Testing surface code...\n";
    stim::CircuitGenParameters surf_code_params(
                            100, 9, "rotated_memory_z");
    surf_code_params.after_clifford_depolarization = 0.003;
    surf_code_params.before_round_data_depolarization = 0.003;
    surf_code_params.before_measure_flip_probability = 0.003;
    surf_code_params.after_reset_flip_probability = 0.003;
    stim::Circuit surf_code_circ =
        generate_surface_code_circuit(surf_code_params).circuit;
    stim::simd_bit_table surf_code_samples = 
        stim::detector_samples(surf_code_circ, 1, false, false, RNG);
    uint8_t surf_code_n_det = 80;  // d^2 - 1 
    for (uint i = 0; i < 100; i++) {
        std::string detection_event;
        for (uint8_t j = 0; j < surf_code_n_det; j++) {
            if (surf_code_samples[i*surf_code_n_det + j][0]) {
                detection_event.push_back('!');
            } else {
                detection_event.push_back('_');
            }
        }
        std::cout << "Detection event: " << detection_event << "\n";
    }
}

void test_decoding_build() {
    std::cout << "======= DECODING (MWPM) BUILD TEST =======\n";
    uint d = 3;
    fp_t p = 0.001;
    uint32_t shots = 100000;

    stim::CircuitGenParameters surf_code_params(
                            d, d, "rotated_memory_z");
    surf_code_params.after_clifford_depolarization = p;
    surf_code_params.before_round_data_depolarization = p;
    surf_code_params.before_measure_flip_probability = p;
    surf_code_params.after_reset_flip_probability = p;
    stim::Circuit surf_code_circ =
        generate_surface_code_circuit(surf_code_params).circuit;
    DecodingGraph graph = to_decoding_graph(surf_code_circ);
    // Print some statistics in the graph.
    std::cout << "number of nodes in decoding graph: " 
        << boost::num_vertices(graph.base) << "\n";
    std::cout << "number of edges in decoding graph: "
        << boost::num_edges(graph.base) << "\n";
    // Time to decode some errors.
    MWPMDecoder decoder(surf_code_circ); 
    b_decoder_ler(&decoder, shots, RNG);
    std::cout << "number of logical errors: " 
        << decoder.n_logical_errors << "\n";
    std::cout << "logical error rate: " 
        << decoder.n_logical_errors/(fp_t)shots << "\n";
    std::cout << "time taken: "
        << min(decoder.execution_times) << " --> " 
        << max(decoder.execution_times) 
        << "(mean = " << mean(decoder.execution_times)
        << ")\n";
    uint32_t n_total_fast;
    auto fast_syndromes = 
        decoder.nonzero_syndromes_completed_within(1000, n_total_fast);
    std::cout << "number of syndromes completed within 1 us: "
        << fast_syndromes.size() << " of " 
        << n_total_fast << "\n";
    std::cout << "max serviceable Hamming weight: "
        << decoder.max_hamming_wgt_completed_within(1000) 
        << "\n";
}

void test_lilliput_build() {
    std::cout << "======= LILLIPUT BUILD TEST =======\n";
    uint d = 3;
    fp_t p = 0.003;
    uint32_t shots = 10000;

    stim::CircuitGenParameters setup_params(
                            2, d, "rotated_memory_z");
    setup_params.after_clifford_depolarization = p;
    setup_params.before_round_data_depolarization = p;
    setup_params.before_measure_flip_probability = p;
    setup_params.after_reset_flip_probability = p;
    stim::Circuit setup_circ =
        generate_surface_code_circuit(setup_params).circuit;

    MWPMDecoder base(setup_circ); 

    LILLIPUTParams params = {
        250e6,                  // Clock frequency
        false,                  // Is DRAM
        (uint64_t)(1 * GB),     // Memory size
        0,                      // tCAS
        0,                      // tRC
        0,                      // tRP
        0,                      // tRCD
        0,                      // n_channels
        0,                      // n_ranks
        0,                      // n_banks
        0,                      // row_size
        0,                      // line_size
        7,                      // sram_access_time
        d*d-1,                  // detectors_per_round
    };
    stim::CircuitGenParameters surf_code_params(
                            d*3, d, "rotated_memory_z");
    surf_code_params.after_clifford_depolarization = p;
    surf_code_params.before_round_data_depolarization = p;
    surf_code_params.before_measure_flip_probability = p;
    surf_code_params.after_reset_flip_probability = p;
    stim::Circuit surf_code_circ =
        generate_surface_code_circuit(surf_code_params).circuit;
    LILLIPUT decoder = LILLIPUT(surf_code_circ, &base, params, RNG);
    b_decoder_ler(&decoder, shots, RNG);
    std::cout << "number of logical errors: " 
        << decoder.n_logical_errors << "\n";
    std::cout << "logical error rate: " 
        << decoder.n_logical_errors/(fp_t)shots << "\n";
    std::cout << "time taken: "
        << min(decoder.execution_times) << " --> " 
        << max(decoder.execution_times) << "\n";
    std::cout << "memory overhead: "
        << mean(decoder.memory_overheads) << "\n";
    std::cout << "LUT misses: "
        << decoder.lut_misses << "\n";
}
