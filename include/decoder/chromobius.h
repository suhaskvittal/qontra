/*
 *  author: Suhas Vittal
 *  date:   23 December 2023
 *
 *  Wrapper for Gidney and Jones's Chromobius decoder.
 * */

#ifndef QONTRA_CHROMOBIUS_h
#define QONTRA_CHROMOBIUS_h

#include "decoder.h"

#include <stim.h>
#include <chromobius/decode/decoder.h>

namespace qontra {

chromobius::Decoder init_chromobius(stim::Circuit circuit) {
    const chromobius::DecoderConfigOptions options;

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
    return chromobius::Decoder::from_dem(dem, options);
}

class Chromobius : public Decoder {
public:
    Chromobius(DetailedStimCircuit circuit)
        :Decoder(circuit, graph::DecodingGraph::Mode::DO_NOT_BUILD),
        backing_decoder(init_chromobius(circuit))
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref syndrome) override {
        const uint det = circuit.count_detectors(),
                    obs = circuit.count_observables();
        timer.clk_start();
        stim::simd_bits corr = 
            backing_decoder.decode_detection_events({syndrome.u8, syndrome.u8 + syndrome.u8.num_bits_padded()});
        fp_t t = (fp_t) timer.clk_end();

        return (Decoder::result_t) {
            t,
            corr,
            is_error(corr, syndrome)
        };
    }

    std::string name(void) override { return "Chromobius"; }
private:
    chromobius::Decoder backing_decoder;
};

}   // qontra

#endif  // QONTRA_CHROMOBIUS_h
