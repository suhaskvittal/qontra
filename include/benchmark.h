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
b_decoder_ler(Decoder*, uint32_t shots, std::mt19937_64&,
        bool save_per_shot_data=true);

}  // qrc

#endif
