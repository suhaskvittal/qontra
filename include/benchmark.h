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

namespace qrc {

fp_t 
min(const std::vector<fp_t>&);
fp_t
max(const std::vector<fp_t>&);
fp_t
mean(const std::vector<fp_t>&);
fp_t
stdev(const std::vector<fp_t>&);

void
b_decoder_ler(Decoder*, uint64_t shots, std::mt19937_64&,
        bool save_per_shot_data=true);


stim::Circuit
build_circuit(
    // Required
    uint code_dist, 
    fp_t error_mean,
    fp_t error_stddev,
    // Optionals
    bool is_memory_z=true,
    bool is_rotated=true,
    bool both_stabilizers=false,
    uint8_t swap_lru=0b00,
    // Level 1 Specificity
    uint rounds=0,
    fp_t clevel_error_mean=-1,
    fp_t clevel_error_stddev=-1,
    fp_t pauliplus_error_mean=-1,
    fp_t pauliplus_error_stddev=-1,
    // Level 2 Specificity
    fp_t round_dp_mean=-1,
    fp_t clifford_dp_mean=-1,
    fp_t reset_flip_mean=-1,
    fp_t meas_flip_mean=-1,
    fp_t round_dp_stddev=-1,
    fp_t clifford_dp_stddev=-1,
    fp_t reset_flip_stddev=-1,
    fp_t meas_flip_stddev=-1,
    fp_t round_leak_mean=-1,
    fp_t clifford_leak_mean=-1,
    fp_t reset_leak_mean=-1,
    fp_t round_leak_stddev=-1,
    fp_t clifford_leak_stddev=-1,
    fp_t reset_leak_stddev=-1);

}  // qrc

#endif
