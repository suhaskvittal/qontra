/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 *
 *  An implementation of Kubica and Delfosse's
 *  Restriction Decoder for color codes.
 * */

#ifndef RESTRICTION_DECODER_h
#define RESTRICTION_DECODER_h

#include "decoder/mwpm.h"
#include "defs.h"
#include "graph/graph.h"

#include <deque>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace qontra {
namespace decoder {

#define __COLOR graph::decoding::vertex_t::Color

inline int c2i(__COLOR x) {
    if (x == __COLOR::red)      return 0;
    if (x == __COLOR::green)    return 1;
    if (x == __COLOR::blue)     return 2;
    std::cout << "NONE OF THE COLORS?\n";
    return -1;
}

class RestrictionDecoder : public Decoder {
public:
    RestrictionDecoder(const stim::Circuit&);

    ~RestrictionDecoder(void) {
        for (auto dec : rlatt_dec)  delete dec;
    }

    std::string name(void) override { return "RestrictionDecoder"; }

    Decoder::result_t   decode_error(const syndrome_t&) override;
private:
    // The decode_restricted_lattice function decodes the restricted lattice that
    // restricts the passed in color.
    Decoder::result_t   decode_restricted_lattice(const std::vector<uint>&, __COLOR);

    typedef std::pair<uint, __COLOR>    cdet_t;
    uint64_t cdet_to_id(cdet_t x) {
        uint64_t base = x.first;
        if (base == graph::BOUNDARY_INDEX) base = circuit.count_detectors();
        return base | ((uint64_t)c2i(x.second) << 60);
    }

    std::array<MWPMDecoder*, 3> rlatt_dec;  // Three MWPM decoders, each of which responsible
                                            // for decoding a different restricted lattice.
                                            //
                                            // Note RestrictionDecoder is a friend of MWPMDecoder,
                                            // so the decoder can access their DecodingGraphs.
    typedef std::map<uint, std::array<uint, 3>>     fwdetmap_t;
    typedef std::map<std::pair<uint, uint>, uint>   bwdetmap_t;
    // As the MWPM decoders operate on
    // different decoding graphs and these
    // graphs do not match the input
    // circuit, we must be able to convert
    // the input syndrome to syndromes for
    // each MWPM decoder.
    fwdetmap_t  to_rlatt; 
    bwdetmap_t  from_rlatt;

    std::map<uint, __COLOR>   color_map;
};

}   // decoder
}   // qontra

#endif  // RESTRICTION_DECODER_h

