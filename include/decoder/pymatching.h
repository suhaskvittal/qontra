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

inline pm::Mwpm init_solver_from_circuit(stim::Circuit circuit) {
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

class PyMatching : public Decoder {
public:
    PyMatching(DetailedStimCircuit circuit)
        :Decoder(circuit, graph::DecodingGraph::Mode::DO_NOT_BUILD),
        solver(init_solver_from_circuit(circuit))
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) override {
        const uint n_observables = circuit.count_observables();

        std::vector<uint> detectors = get_nonzero_detectors_(syndrome);
        if (detectors.size() == 0) {
            stim::simd_bits<SIMD_WIDTH> corr(n_observables);
            corr.clear();
            return (Decoder::result_t) {
                0,
                corr
            };
        }
        // Convert to uint64_t.
        std::vector<uint64_t> d64;
        for (uint x : detectors) d64.push_back(x);
        stim::simd_bits<SIMD_WIDTH> corr(n_observables);

        // Get working set into the cache first.
        int64_t w;
        pm::decode_detection_events(solver, d64, corr.u8, w);
        corr.clear();
        // Now time for realz.
        timer.clk_start();
        pm::decode_detection_events(solver, d64, corr.u8, w);
        fp_t t = (fp_t)timer.clk_end();

        return (Decoder::result_t) {
            t,
            corr
        };
    }

    inline std::string name(void) override { return "PyMatching"; }
private:
    pm::Mwpm solver;
};

}   // qontra

#endif  // QONTRA_PYMATCHING_h
