/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef DECODER_h
#define DECODER_h

#include <stim.h>

#include "defs.h"
#include "defs/timer.h"
#include "graph/decoding_graph.h"
#include "stimext.h"

#include <string>
#include <tuple>

#include <strings.h>

namespace qontra {

template <class T> std::vector<uint>
// T = stim::simd_bits or stim::simd_bits_range_ref
get_nonzero_detectors(T syndrome, uint number_of_detectors) {
    std::vector<uint> det;
    uint64_t w = 0;
    uint64_t last_bit = 0;
    while (det.size() < syndrome.popcnt()) {
        uint64_t i = ffsll(syndrome.u64[w] & ~((1L << last_bit)-1));
        if (i == 0) {   // No match found.
            last_bit = 0;
            w++;
            continue;
        }
        uint d = (w << 6) | (i-1);
        if (d >= number_of_detectors) break;
        det.push_back(d);
        last_bit = i & 0x3f;
        w += (i >= 64);
    }
    return det;
}

class Decoder {
public: 
    // The decoder uses the passed in circuit for decoding. It should
    // be a memory experiment circuit for X rounds.
    //
    // The second argument in the constructor is optional and prevents
    // the decoding graph from being constructed. The user should not
    // use this, but rather a subclass constructor should.
    // This is useful for sliding window decoders, where
    // the circuit used for evaluation may be 100s or 1000s of rounds long,
    // but the backing decoder only uses a few rounds. This may also be 
    // used if the subclass constructs the decoding graph
    // via a custom strategy, or does not need a decoding graph.
    Decoder(DetailedStimCircuit circ, 
            graph::DecodingGraph::Mode dgr_mode=graph::DecodingGraph::Mode::NORMAL)
        :circuit(circ),
        decoding_graph(graph::to_decoding_graph(circ, dgr_mode))
    {}

    virtual ~Decoder() {}

    typedef std::tuple<uint, uint, stim::simd_bits<SIMD_WIDTH>> assign_t;

    struct result_t {
        fp_t exec_time = 0.0;
        stim::simd_bits<SIMD_WIDTH> corr = stim::simd_bits<>(1);
        // The below arguments are not guaranteed to be populated.
        std::vector<assign_t>   error_assignments;
    };

    virtual result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) =0;
                                                    // This function
                                                    // will perform decoding
                                                    // and return a correction.
    virtual std::string name(void) =0;  // Useful for printing out stats.

    DetailedStimCircuit get_circuit(void);
protected:
    template <class T> std::vector<uint> get_nonzero_detectors_(T syndrome);

    DetailedStimCircuit     circuit;
    graph::DecodingGraph    decoding_graph;

    Timer   timer;

    friend class WindowDecoder;
};

}   // qontra

#include "decoder.inl"

#endif
