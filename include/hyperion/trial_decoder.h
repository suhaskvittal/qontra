/*
 *  author: Suhas Vittal
 *  date:   1 October 2022
 * */

#ifndef HYPERION_TRIAL_DECODER_h
#define HYPERION_TRIAL_DECODER_h

#include "defs.h"
#include "mwpm_decoder.h"

#include <deque>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

#include <time.h>

#define BDC_DEBUG

namespace qrc {
namespace hyperion {

class Blossom {
public:
    Blossom(uint node_id);

    bool contains(Blossom*);
    void reset(void);

    enum class Label {S, T, Free};

    uint node_id;
    bool is_blossom; 
    Label label;
    std::vector<Blossom*> bfwdptrs;
    Blossom * bbwdptr;
    bool is_visited;
    bool is_matched;
    Blossom * mateptr;
    Blossom * prevptr;
    Blossom * ownerptr;

    wgt_t u;
    wgt_t z;
}; 

class TrialDecoder : public MWPMDecoder {
public:
    TrialDecoder(const stim::Circuit&);
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
private:
    std::map<uint, uint> hungarian(const std::vector<uint>&);
    std::map<uint, uint> blossom(const std::vector<uint>&);
};

}   // gulliver
}   // qrc

#endif
