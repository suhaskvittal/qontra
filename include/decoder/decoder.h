/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef DECODER_h
#define DECODER_h

#include <stim.h>

#include "defs.h"
#include "graph/decoding_graph.h"

#include <string>
#include <utility>

namespace qontra {
namespace decoder {

class Decoder {
public: 
    Decoder(const stim::Circuit& circ)          // The decoder uses the passed in
                                                // circuit for decoding. It should
                                                // be a memory experiment circuit
                                                // for X rounds.
        :circuit(circ),
        decoding_graph(graph::to_decoding_graph(circ))
    {}

    virtual ~Decoder() {}

    typedef struct {
        fp_t exec_time;
        std::vector<uint8_t> corr;
        bool is_error;
    } result_t;

    typedef std::vector<uint8_t>    vsyndrome_t;

    virtual result_t decode_error(const vsyndrome_t&) =0;   // This function
                                                            // will perform decoding
                                                            // and return a correction.
    virtual std::string name(void) =0;  // Useful for printing out stats.

    stim::Circuit       get_circuit(void) { return circuit; }
protected:
    // Other helpful functions:
    //
    // clk_start: records the time it was called.
    // clk_end: records the time elapsed since the last clk_start.
    // is_error: checks if the provided correction is a logical error.

    void clk_start(void) {
#ifdef __APPLE__
        t_start = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
        struct timespec d;
        clock_gettime(CLOCK_MONOTONIC_RAW, &d);
        t_start = d.tv_nsec;
#endif
    }

    uint64_t clk_end(void) {
        auto tmp = t_start;
        clk_start();
        return t_start - tmp;
    }

    bool is_error(
            const std::vector<uint8_t>& correction, const vsyndrome_t& syndrome)
    {
        const uint n_detectors = circuit.count_detectors();
        const uint n_observables = circuit.count_observables();
        bool is_error = false;
        for (uint i = 0; i < n_observables; i++) {
            is_error |= correction[i] ^ syndrome[n_detectors+i];
        }
        return is_error;
    }

    stim::Circuit circuit;
    graph::DecodingGraph decoding_graph;
private:
#ifdef __APPLE__
    uint64_t    t_start;
#else
    long        t_start;
#endif
};

// Helper functions
// 
// syndrome_to_vector converts a simd_bits_range_ref to a vector. A vector
// can be modified and is easier to pass around.

inline Decoder::vsyndrome_t
syndrome_to_vector(const stim::simd_bits_range_ref& ref, uint size) {
    std::vector<uint8_t> syndrome(size);
    for (uint i = 0; i < size; i++) syndrome[i] = ref[i];
    return syndrome;
}

}   // decoder
}   // qontra

#endif
