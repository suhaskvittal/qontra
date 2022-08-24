/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef GULLIVER_EXPERIMENTS_h
#define GULLIVER_EXPERIMENTS_h

#include "quarch.h"
#include "quarch_gulliver.h"
#include "quarch_defs.h"

#include <chrono>
#include <filesystem>
#include <random>
#include <string>

extern std::filesystem::path data_folder;
extern std::mt19937_64 GULLIVER_RNG;

struct GulliverTimingAnalysisParams {
    Decoder * decoder_p;
    uint32_t shots;
};

struct GulliverSweepAnalysisParams {
    ErrorThresholdSweepParams::decoder_gen_f dgf;

    fp_t error_lb;
    fp_t error_ub;
    fp_t error_step;

    uint code_distance_lb;
    uint code_distance_ub;
    uint code_distance_step;

    uint32_t shots;
};

void 
decoder_analysis_experiment();
void
gulliver_timing_experiment();
void 
mwpm_timing_experiment();
void
mwpm_sweep_experiment();

/* Private Helper Functions */
void
_timing_analysis(const std::filesystem::path& output_file, 
        const GulliverTimingAnalysisParams&);
void
_sweep_analysis(const std::filesystem::path& output_folder,
        const GulliverSweepAnalysisParams&);

#endif
