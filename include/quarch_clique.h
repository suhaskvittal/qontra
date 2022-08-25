/*
 *  author: Suhas Vittal
 *  date:   25 August 2022
 * */

#ifndef QUARCH_CLIQUE_h
#define QUARCH_CLIQUE_h

#include "quarch_mwpm_decoder.h"

class CliqueDecoder : public MWPMDecoder: {
public:
    CliqueDecoder(const stim::Circuit&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&); override;
    std::string name(void) override;
    bool is_software(void) override;

    // More statistics.
    uint32_t complex_decoder_uses;
};

#endif
