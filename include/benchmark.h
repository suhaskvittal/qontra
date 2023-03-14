/*
 *  author: Suhas Vittal
 *  date:   21 August 2022
 * */

#ifndef BENCHMARK_h
#define BENCHMARK_h

#include <stim.h>

#include "benchmark/statbench.h"
#include "benchmark/statbench/analytical_dist.h"
#include "benchmark/statbench/numerical_dist.h"
#include "defs.h"
#include "decoder.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>

#include <math.h>
#include <mpi.h>

#define MAX_SHOTS 100000

namespace qrc {

fp_t 
min(const std::vector<fp_t>&);
fp_t
max(const std::vector<fp_t>&);
fp_t
mean(const std::vector<fp_t>&);
fp_t
stddev(const std::vector<fp_t>&);

void
b_decoder_ler(Decoder*, uint64_t shots, std::mt19937_64&, bool save_per_shot_data=false);
/*
 *  Pre-condition: MPI is initialized before call and exited after call.
 * */
benchmark::StatisticalResult
b_statistical_ler(
    dgf_t&,
    uint code_dist,
    fp_t start_p,
    fp_t final_p, 
    uint64_t shots,
    uint64_t update_rate,
    std::mt19937_64&, 
    bool use_mpi=false,
    bool bootstrap_model=false, 
    std::map<uint, uint64_t> bootstrap_data=std::map<uint, uint64_t>(),
    fp_t use_bootstrap_model_until_p=1
);

#define BC_FLAG_SWAP_LRU_V1     0x1
#define BC_FLAG_SWAP_LRU_V2     0x2
#define BC_FLAG_INVERT_STATE    0x4

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
    uint8_t other_flags=0,
    // Level 1 Specificity
    uint rounds=0,
    fp_t clevel_error_mean=-1,
    fp_t clevel_error_stddev=-1,
    fp_t pauliplus_error_mean=-1,
    fp_t pauliplus_error_stddev=-1,
    // Level 2 Specificity
    fp_t round_dp_mean=-1,
    fp_t sq_dp_mean=-1,
    fp_t cx_dp_mean=-1,
    fp_t reset_flip_mean=-1,
    fp_t meas_flip_mean=-1,
    fp_t round_dp_stddev=-1,
    fp_t sq_dp_stddev=-1,
    fp_t cx_dp_stddev=-1,
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
