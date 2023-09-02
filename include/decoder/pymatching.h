/*
 *  author: Suhas Vittal
 *  date:   2 September 2023
 *
 *  A simple compatability wrapper
 *  for PyMatching.
 * */

#ifndef QONTRA_PYMATCHING_h
#define QONTRA_PYMATCHING_h

#include "decoder/decoder.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

namespace qontra {

pm::Mwpm init_solver_from_circuit(stim::Circuit circuit) {
    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit,
            true,  // decompose_errors
            true,  // fold loops
            true, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            true, // ignore decomposition failures
            false
        );
    return pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS);
}

class PyMatching : public Decoder {
public:
    PyMatching(stim::Circuit circuit)
        :Decoder(circuit, graph::DecodingGraph::Mode::LOW_MEMORY),
        solver(init_solver_from_circuit(circuit))
    {}

    Decoder::result_t decode_error(const syndrome_t& syndrome) override {
        std::vector<uint> detectors = get_nonzero_detectors(syndrome);
        // Convert to uint64_t.
        std::vector<uint64_t> d64;
        for (uint x : detectors) d64.push_back(x);
        
        const uint n_observables = circuit.count_observables();
        stim::simd_bits corr(n_observables);

        timer.clk_start();
        int64_t w;
        pm::decode_detection_events(solver, d64, corr.u8, w);
        fp_t t = (fp_t)timer.clk_end();

        return (Decoder::result_t) {
            t,
            corr,
            is_error(corr, syndrome)
        };
    }

    std::string name(void) override { return "PyMatching"; }
private:
    pm::Mwpm solver;
};

}   // qontra

#endif  // QONTRA_PYMATCHING_h
