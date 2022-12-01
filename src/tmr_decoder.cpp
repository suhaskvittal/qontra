/*
 *  author: Suhas Vittal
 *  date:   1 December 2022
 * */

#include "tmr_decoder.h"

namespace qrc {

TMRDecoder::TMRDecoder(
        const stim::Circuit& circuit,
        Decoder * baseline,
        uint detectors_per_round)
    :Decoder(circuit),
    baseline(baseline),
    detectors_per_round(detectors_per_round)
{}

std::string
TMRDecoder::name() {
    return (baseline->name() + "[TMR]");
}

bool
TMRDecoder::is_software() {
    return baseline->is_software();
}

DecoderShotResult
TMRDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    uint tmr_detectors = n_detectors / 3;
    std::vector<uint8_t> tmr_syndrome(tmr_detectors + n_observables);

    std::vector<uint8_t> ssyndrome(syndrome);
    unxor(ssyndrome);
    // Convert unxor'd syndrome to TMR syndrome.
    uint8_t last_bit = 0;
    for (uint i = 0; i < tmr_detectors; i++) {
        uint k1 = 3*i, k2 = 3*i + 1, k3 = 3*i + 2;
        uint8_t bit = (ssyndrome[k1] + ssyndrome[k2] + ssyndrome[k3]) >= 2;
        tmr_syndrome[i] = bit ^ last_bit;
        last_bit = bit;
    }
    // Also add observables (unchanged).
    for (uint i = 0; i < n_observables; i++) {
        tmr_syndrome[tmr_detectors + i] = ssyndrome[n_detectors + i];
    }

    return baseline->decode_error(tmr_syndrome);
}

void
TMRDecoder::unxor(std::vector<uint8_t>& syndrome) {
    uint n_observables = circuit.count_observables();
    for (uint i = 1; i < syndrome.size() - n_observables; i++) {
        syndrome[i] ^= syndrome[i-1];
    }
}

}   // qrc
