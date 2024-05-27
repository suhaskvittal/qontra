/*
 *  author: Suhas Vittal
 *  date:   22 January 2024
 * */

#include "qontra/decoder/pymatching.h"

namespace qontra {

pm::Mwpm
init_solver_from_circuit(stim::Circuit circuit) {
    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit,
            true,  // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    return pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS);
}

PyMatching::PyMatching(const DetailedStimCircuit& circuit)
    :Decoder(circuit),
    solver(init_solver_from_circuit(circuit))
{}

Decoder::result_t
PyMatching::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_observables = circuit.count_observables();

    std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
    stim::simd_bits<SIMD_WIDTH> corr(n_observables);

    int64_t w;
    timer.clk_start();
    pm::decode_detection_events(solver, detectors, corr.u8, w);
    fp_t t = (fp_t)timer.clk_end();

    return { t, corr };
}

}   // qontra

