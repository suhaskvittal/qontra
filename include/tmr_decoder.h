/*
 *  author: Suhas Vittal
 *  date:   1 December 2022
 * */

#ifndef TMR_DECODER_h
#define TMR_DECODER_h

#include "decoder.h"
#include "defs.h"

namespace qrc {

class TMRDecoder : public Decoder {
public:
    TMRDecoder(const stim::Circuit&, Decoder*, uint detectors_per_round);

    std::string name(void) override;
    bool is_software(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
private:
    void unxor(std::vector<uint8_t>&);

    Decoder * baseline;
    uint detectors_per_round;
};


}   // qrc

#endif
