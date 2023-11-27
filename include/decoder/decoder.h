/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef DECODER_h
#define DECODER_h

#include <stim.h>

#include "defs.h"
#include "graph/decoding_graph.h"
#include "stimext.h"
#include "timer.h"

#include <string>
#include <tuple>
#include <utility>

#include <strings.h>

namespace qontra {

template <class T> std::vector<uint>
// T = stim::simd_bits or stim::simd_bits_range_ref
get_nonzero_detectors(const T syndrome, uint number_of_detectors) {
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

    typedef std::tuple<uint, uint, stim::simd_bits> assign_t;

    struct result_t {
        fp_t exec_time = 0.0;
        stim::simd_bits corr = stim::simd_bits(1);
        bool is_error = false;

        std::vector<assign_t>   error_assignments;
    };

    virtual result_t decode_error(stim::simd_bits_range_ref) =0;
                                                    // This function
                                                    // will perform decoding
                                                    // and return a correction.
    virtual std::string name(void) =0;  // Useful for printing out stats.

    DetailedStimCircuit get_circuit(void) { return circuit; }
protected:
    // Other helpful functions:
    //
    // is_error: checks if the provided correction is a logical error.
    // get_nonzero_detectors: gets nonzero detectors from syndrome.

    bool is_error(
            const stim::simd_bits& correction, stim::simd_bits_range_ref syndrome)
    {
        const uint n_detectors = circuit.count_detectors();
        const uint n_observables = circuit.count_observables();
        bool is_error = false;
        for (uint i = 0; i < n_observables; i++) {
            is_error |= correction[i] ^ syndrome[n_detectors+i];
        }
        return is_error;
    }

    template <class T> std::vector<uint>
    // T = stim::simd_bits or stim::simd_bits_range_ref
    get_nonzero_detectors(T syndrome) {
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
            if (d >= circuit.count_detectors()) break;
            det.push_back(d);
            last_bit = i & 0x3f;
            w += (i >= 64);
        }
        return det;
    }

    DetailedStimCircuit     circuit;
    graph::DecodingGraph    decoding_graph;

    Timer   timer;

    friend class WindowDecoder;
};

}   // qontra
#endif
