/*
 *  author: Suhas Vittal
 *  date:   21 August 2022
 * */

#ifndef BENCHMARK_h
#define BENCHMARK_h

#include <stim.h>

#include "defs.h"
#include "decoder.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <utility>

#define MAX_SHOTS 100000

struct ErrorThresholdData {
    // Key: (code distance, error rate)
    // Value: (circuit, logical error rate)
    typedef std::pair<uint, fp_t> DataKey;
    typedef std::pair<stim::Circuit, fp_t> DataValue;

    std::map<DataKey, DataValue> data;
};

struct ErrorThresholdSweepParams {
    typedef std::function<Decoder*(uint, fp_t)> decoder_gen_f;
    // Function will create new decoder for code distance,
    // physical error rate pair.
    // We will assume that the output pointer points to the heap.
    decoder_gen_f decoder_generator;
    // All upper bounds are inclusive.
    fp_t error_lb;
    fp_t error_ub;
    fp_t error_step;

    uint code_distance_lb;
    uint code_distance_ub;
    uint code_distance_step;
};

struct ErrorThresholdOutputParams {
    bool save_circuits_to_folder;
    std::string circuit_folder_path;
};

fp_t 
min(const std::vector<fp_t>&);
fp_t
max(const std::vector<fp_t>&);
fp_t
mean(const std::vector<fp_t>&);
fp_t
stdev(const std::vector<fp_t>&);

void
b_decoder_ler(Decoder*, uint32_t shots, std::mt19937_64&,
        bool save_per_shot_data=true);

ErrorThresholdData
sweep_error_threshold(const ErrorThresholdSweepParams&,
        uint32_t shots, std::mt19937_64& rng);
void
write_sweep_data(const ErrorThresholdData&, std::ostream&,
        const ErrorThresholdOutputParams&);

#endif
