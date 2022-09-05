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

/*
 * Determine frequency of various Hamming
 * weights for different distances
 * (surface code).
 * */
void
surface_code_hamming_weight_experiment();

/* Private Helper Functions */
stim::Circuit
_make_surface_code_circuit(uint code_dist, fp_t, fp_t);
void
_timing_analysis(const std::filesystem::path& output_file, 
        Decoder*, uint32_t shots);

#endif
