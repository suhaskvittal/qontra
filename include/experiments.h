/* author: Suhas Vittal
 *  date:   12 June 2023
 * */

#ifndef EXPERIMENTS_h
#define EXPERIMENTS_h

#include "decoder/decoder.h"
#include "defs.h"

#include <stim.h>

#include <functional>
#include <string>

#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

namespace qontra {

namespace experiments {

// We define two types of callbacks:
//
// cb_t1 is a callback that operates on the syndrome (simd_bits_range_ref). cb_t1
// is a generic callback that can be used with generate_syndromes (see below).
//
// cb_t2 is a callback that operates on a decoder's result. cb_t2 works with
// memory experiments and should be used to update stats.
//
// Both callbacks are packaged into the callback_t struct.

typedef std::function<void(stim::simd_bits_range_ref&)>         cb_t1;
typedef std::function<void(const Decoder::result_t&)>  cb_t2;

typedef struct {
    cb_t1   prologue = [] (stim::simd_bits_range_ref x) {};
    cb_t2   epilogue = [] (const Decoder::result_t& res) {};
} callback_t;

extern callback_t   DEFAULT_CALLBACKS;

// The memory_params_t struct contains all parameters to execute a memory
// experiment. The user can choose to use their own trace, or to run
// the experiment on generated syndromes.
//  (1) If using a trace, then set shots = 0.
//  (2) Else, set trace_folder = "DNE".

typedef struct {
    uint64_t    shots           = 0;
    std::string trace_folder    = "DNE";

    callback_t  callbacks = DEFAULT_CALLBACKS;
} memory_params_t;

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

// Useful defines:

#define SQR(x)              (x)*(x)
#define MEAN(s, n)          ((fp_t)(s))/((fp_t)(n))
#define STD(m, ss, n)       sqrt( ((fp_t)(ss))/((fp_t)(n)) - SQR(m) )

// Experiments functions:
//
// generate_syndromes generates syndromes via Monte carlo sampling and
// executes the given callback for each syndrome.
//
// build_syndrome_trace generates syndromes and writes them to a folder
// in the .dets format (from Stim).
//
// memory_experiment calls generate_syndromes by supplying a callback that
// calls the provided decoder. The function is overloaded to also provide
// a cb_t1 parameter that can be used in conjunction with the generated cb_t1
// callback for the decoder.
//
void    generate_syndromes(const stim::Circuit&,
                                uint64_t shots,
                                experiments::callback_t);
void    build_syndrome_trace(std::string, const stim::Circuit&, uint64_t shots);

uint64_t
read_syndrome_trace(std::string, const stim::Circuit&, experiments::callback_t);

experiments::memory_result_t
memory_experiment(Decoder*, experiments::memory_params_t);

}   // qontra

#endif  // EXPERIMENTS_h
