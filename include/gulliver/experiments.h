/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef GULLIVER_EXPERIMENTS_h
#define GULLIVER_EXPERIMENTS_h

#include "quarch.h"
#include "gulliver.h"
#include "clique.h"
#include "defs.h"

#include <chrono>
#include <filesystem>
#include <random>
#include <string>

extern std::filesystem::path data_folder;
extern std::mt19937_64 GULLIVER_RNG;

struct TimingAnalysisParams {
    Decoder * decoder_p;
    uint32_t shots;
};

/* Counts amount of SRAM storage
 * required to support decoders. */
void
decoder_sram_experiment();

/* Examines the logical error rate
 * and syndrome completion of 
 * decoders.
 * */
void 
decoder_analysis_experiment();

/* Computes the timing of various
 * sampled syndromes using various
 * decoers.
 * */
void
gulliver_timing_experiment();
void 
mwpm_timing_experiment();

/* Performs a sweep to get data
 * to compute the error threshold.
 * */
void
mwpm_sweep_experiment();

/* Examines the performance of 
 * Gulliver's cache with different
 * parameters across different
 * distances.
 * */
void
gulliver_cache_experiment();

/* Private Helper Functions */
void
_timing_analysis(const std::filesystem::path& output_file, 
        Decoder*, uint32_t shots);
void
_threshold_sweep(const std::filesystem::path& output_folder,
        const ErrorThresholdSweepParams&, uint32_t shots);

void
_cache_sweep(const std::filesystem::path& output_file,
        const stim::Circuit&, uint32_t shots);

#endif
