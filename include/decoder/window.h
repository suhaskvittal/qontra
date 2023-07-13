/*
 *  author: Suhas Vittal
 *  date:   7 July 2023
 * */

#ifndef DECODER_WINDOW_h
#define DECODER_WINDOW_h

#include "decoder/decoder.h"
#include "decoder/mwpm.h"

namespace qontra {
namespace decoder {

// This is a class for sliding window decoders.
// The sliding window decoder uses a "base" decoder.
// This decoder should be capable of decoding
// 1 round than wanted (d+1, for example). This
// is to avoid time boundary effects.
//
// Corrections are computed within the "commit"
// window. It is assumed that the experiment
// handles how fast the window moves.

class WindowDecoder : public Decoder {
public:
    WindowDecoder(const stim::Circuit& backing_circuit,
                    Decoder* base,
                    uint commit_window, 
                    uint detectors_per_round)
        :Decoder(backing_circuit),
        base_decoder(base),
        commit_window(commit_window),
        detectors_per_round(detectors_per_round)
    {}

    Decoder::result_t   decode_error(const syndrome_t&) override;
                                // This function (which overrides from
                                // the superclass) decodes the entire
                                // provided syndrome over a sliding window.
                                // The syndrome should be for the backing
                                // circuit.
                                //
                                // This is really only for testing.

    // The below functions handle syndromes at different points. The
    // decode_middle_error should be called most often.
    Decoder::result_t   decode_first_error(const syndrome_t&);
    Decoder::result_t   decode_middle_error(const syndrome_t&);
    Decoder::result_t   decode_final_error(const syndrome_t&);

    std::string name(void) override { return "W" + base_decoder->name(); }
protected:
    // This helper function calls the base decoder and retrieves the
    // result such that only error assignments made involving detectors
    // between the min (inclusive) and max (exclusive) are kept.
    Decoder::result_t   retrieve_result_from_base(
                                        const syndrome_t& base_syndrome,
                                        uint min_detector,
                                        uint max_detector);
private:
    Decoder*    base_decoder;
    uint        commit_window;
    uint        detectors_per_round;
};

}   // decoder
}   // qontra

#endif  // DECODER_WINDOW_h
