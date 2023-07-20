/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 *
 *  An implementation of Kubica and Delfosse's
 *  Restriction Decoder for color codes
 *  (with Chamberland et al.'s extension).
 * */

#ifndef RESTRICTION_DECODER_h
#define RESTRICTION_DECODER_h

#include "decoder/mwpm.h"
#include "defs.h"

namespace qontra {
namespace decoder {

// To distinguish the color of a detector, we will
// assume by default that det % 0 =
//      0   if the detector is a "red" event
//      1   if it is blue
//      2   if it is green
// Naturally, the boundary is exempt from this
// requirement.
//
// The user can also specify the color of detector.

#define CC_RED      0
#define CC_BLUE     1
#define CC_GREEN    2

class RestrictionDecoder : public MWPMDecoder {
public:
    RestrictionDecoder(const stim::Circuit& circuit)
        :MWPMDecoder(circuit)
    {
        for (uint64_t d = 0; d < circuit.count_detectors(); d++) {
            auto v = decoding_graph.get_vertex(d);
            if (d % 3 == CC_RED)        detector_to_color[v] = Color::red;
            else if (d % 3 == CC_BLUE)  detector_to_color[v] = Color::blue;
            else if (d % 3 == CC_GREEN) detector_to_color[v] = Color::green;
        }
    }

    enum class Color { red, blue, green };

    std::string name(void) override { return "RestrictionDecoder"; }

    void set_detector_color(uint d, Color c) {
        auto v = decoding_graph.get_vertex(d);
        detector_to_color[v] = c;
    }

    Decoder::result_t   decode_error(const syndrome_t&) override;
private:
    Decoder::result_t   decode_restricted_lattice(
                            const std::vector<uint>& detectors,
                            Color c1,
                            Color c2);

    std::map<graph::decoding::vertex_t*, Color> detector_to_color;
};

}   // decoder
}   // qontra

#endif  // RESTRICTION_DECODER_h

