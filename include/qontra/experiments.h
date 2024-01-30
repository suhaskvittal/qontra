/* 
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#ifndef QONTRA_EXPERIMENTS_h
#define QONTRA_EXPERIMENTS_h

#include "qontra/ext/stim.h"

#include <string>

namespace qontra {

struct shot_payload_t {
    stim::simd_bits_range_ref<SIMD_WIDTH>   syndrome;
    stim::simd_bits_range_ref<SIMD_WIDTH>   observables;
};

// Global experiments parameters:
extern bool     G_USE_MPI;                          // Default is true.
extern uint64_t G_SHOTS_PER_BATCH;                  // Default is 100'000.
extern uint64_t G_BASE_SEED;                        // Default is 0.
extern bool     G_FILTER_OUT_SYNDROMES;             // Default is true.
extern uint64_t G_FILTERING_HAMMING_WEIGHT;         // Default is 0.

void    configure_optimal_batch_size(void);

std::string get_batch_filename(size_t batchno);

template <class CALLBACK>
void generate_syndromes(const DetailedStimCircuit&, uint64_t shots, CALLBACK);

template <class CALLBACK>
uint64_t read_syndrome_trace(std::string input_folder, const DetailedStimCircuit&, CALLBACK);

void build_syndrome_trace(std::string output_folder, const DetailedStimCircuit&, uint64_t shots);

}   // qontra

#include "experiments.inl"

#endif  // QONTRA_EXPERIMENTS_h
