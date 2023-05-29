/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef DECODER_h
#define DECODER_h

#include <stim.h>

#include "defs.h"
#include "decoding_graph.h"

#include <string>
#include <utility>

namespace qontra {

class Decoder {
public: 
    Decoder(const stim::Circuit& circ)          // The decoder uses the passed in
                                                // circuit for decoding. It should
                                                // be a memory experiment circuit
                                                // for X rounds.
        :circuit(circ),
        decoding_graph(to_decoding_graph(circ))
    {}

    virtual ~Decoder() {}

    typedef struct {
        fp_t exec_time;
        std::vector<uint8_t> corr;
    } result_t;

    typedef std::vector<uint8_t>    vsyndrome_t;

    virtual result_t decode_error(const vsyndrome_t&) =0;   // This function
                                                            // will perform decoding
                                                            // and return a correction.
    virtual std::string name(void) =0;  // Useful for printing out stats.
protected:
    // Other helpful functions:
    void clk_start(void) {  // Records the time the function was called.
#ifdef __APPLE__
        t_start = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
        struct timespec d;
        clock_gettime(CLOCK_MONOTONIC_RAW, &d);
        t_start = d.tv_nsec;
#endif
    }

    uint64_t clk_end(void) {    // Returns the time elapsed since the last clk_start
                                // call.
        auto tmp = t_start;
        clk_start();
        return tmp - t_start;
    }

    stim::Circuit circuit;
    DecodingGraph decoding_graph;
private:
#ifdef __APPLE__
    uint64_t    t_start;
#else
    long        t_start;
#endif
};

// Helper functions
// 
// to_vector converts a simd_bits_range_ref to a vector. A vector
// can be modified and is easier to pass around.

std::vector<uint8_t> 
to_vector(const stim::simd_bits_range_ref& ref, uint size) {
    std::vector<uint8_t> syndrome(size);
    for (uint i = 0; i < size; i++) syndrome[i] = ref[i];
    return syndrome;
}

}   // qontra

#endif
