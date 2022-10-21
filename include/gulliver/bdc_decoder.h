/*
 *  author: Suhas Vittal
 *  date:   1 October 2022
 * */

#ifndef GULLIVER_BDC_DECODER_h
#define GULLIVER_BDC_DECODER_h

#include "defs.h"
#include "mwpm_decoder.h"

#include <deque>
#include <limits>
#include <map>
#include <vector>

#include <time.h>

#define BDC_DEBUG

namespace qrc {
namespace gulliver {

#define SBOY    0
#define TGIRL   1
#define FREE    2

class BDCDecoder : public MWPMDecoder {
public:
    BDCDecoder(const stim::Circuit&);
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
private:
    std::map<uint, uint> bdc_mwpm(const std::vector<uint>&);
};

}   // gulliver
}   // qrc

#endif
