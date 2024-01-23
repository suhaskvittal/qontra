/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef QONTRA_DECODER_h
#define QONTRA_DECODER_h

#include "qontra/defs.h"
#include "qontra/graph/decoding_graph.h"
#include "qontra/ext/stim.h"

#include <vtils/timer.h>

#include <string>
#include <tuple>

#include <strings.h>

namespace qontra {

template <class T>
// T = stim::simd_bits or stim::simd_bits_range_ref
std::vector<uint64_t> get_nonzero_detectors_(T, uint64_t number_of_detectors);

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

    typedef std::tuple<uint64_t, uint64_t, stim::simd_bits<SIMD_WIDTH>> assign_t;

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

    DetailedStimCircuit get_circuit(void);
protected:
    template <class T>
    std::vector<uint64_t> get_nonzero_detectors(T syndrome);

    DetailedStimCircuit     circuit;
    graph::DecodingGraph    decoding_graph;

    vtils::Timer    timer;

    friend class WindowDecoder;
};

}   // qontra

#include "decoder.inl"

#endif
