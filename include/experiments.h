/*
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#include "decoder/decoder.h"
#include "defs.h"

#include <stim.h>

#include <functional>

#include <mpi.h>

namespace qontra {

namespace experiments {

// We define two types of callbacks:
//
// cb_t1 is a callback that operates on the syndrome (simd_bits_range_ref). cb_t1
// is a generic callback that can be used with generate_syndromes (see below).
//
// cb_t2 is a callback that operates on a decoder's result. cb_t2 works with
// memory experiments and should be used to update stats.

typedef std::function<void(stim::simd_bits_range_ref&)>         cb_t1;
typedef std::function<void(const decoder::Decoder::result_t&)>  cb_t2;

typedef struct {
    fp_t        logical_error_rate;
    fp_t        hw_mean;
    fp_t        hw_std;
    uint64_t    hw_max;
    fp_t        t_mean;
    fp_t        t_std;
    fp_t        t_max;
} memory_result_t;

// Global experiments parameters:
extern bool     G_USE_MPI;                          // Default is true.
extern uint64_t G_SHOTS_PER_BATCH;                  // Default is 100'000.
extern uint64_t G_BASE_SEED;                        // Default is 0.
extern bool     G_FILTER_OUT_SYNDROMES;             // Default is true.
extern uint64_t G_FILTERING_HAMMING_WEIGHT;         // Default is 2.

}   // experiments

// Experiments functions:
//
// generate_syndromes generates syndromes via Monte carlo sampling and
// executes the given callback for each syndrome.
//
// memory_experiment calls generate_syndromes by supplying a callback that
// calls the provided decoder. The function is overloaded to also provide
// a cb_t1 parameter that can be used in conjunction with the generated cb_t1
// callback for the decoder.
//
void    generate_syndromes(const stim::Circuit&, uint64_t shots, experiments::cb_t1);

experiments::memory_result_t
memory_experiment(decoder::Decoder*, uint64_t shots,
                            experiments::cb_t1, experiments::cb_t2);

inline experiments::memory_result_t
memory_experiment(decoder::Decoder* dec, uint64_t shots, experiments::cb_t2 cb) {
    experiments::cb_t1 x = [&] (stim::simd_bits_range_ref&) {};
    return memory_experiment(dec, shots, x, cb);
}

inline experiments::memory_result_t
memory_experiment(decoder::Decoder* dec, uint64_t shots) {
    experiments::cb_t2 x = [&] (decoder::Decoder::result_t res) {};
    return memory_experiment(dec, shots, x);
}

}   // qontra