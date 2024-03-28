/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef QONTRA_DECODER_h
#define QONTRA_DECODER_h

#include "qontra/defs.h"
#include "qontra/ext/stim.h"

#include <vtils/timer.h>

#include <string>
#include <tuple>

namespace qontra {

template <class T>
// T = stim::simd_bits or stim::simd_bits_range_ref
std::vector<uint64_t> get_nonzero_detectors_(T, uint64_t number_of_detectors);

class Decoder {
public: 
    struct result_t {
        fp_t exec_time = 0.0;
        stim::simd_bits<SIMD_WIDTH> corr = stim::simd_bits<>(1);
    };

    // The decoder uses the passed in circuit for decoding. It should
    // be a memory experiment circuit for X rounds.
    Decoder(const DetailedStimCircuit& circ)
        :circuit(circ)
    {}

    virtual result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) =0;
    DetailedStimCircuit get_circuit(void);
    void set_circuit(const DetailedStimCircuit&);
protected:
    template <class T>
    std::vector<uint64_t> get_nonzero_detectors(T syndrome);

    DetailedStimCircuit circuit;
    vtils::Timer timer;
};

}   // qontra

#include "decoder.inl"

#endif
