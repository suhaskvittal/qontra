/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QUARCH_STIM_TEST_h
#define QUARCH_STIM_TEST_h

#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>

#include <stdint.h>

#include <stim.h>

#include "quarch_defs.h"
#include "quarch_decoding_graph.h"
#include "quarch_mwpm_decoder.h"
#include "quarch_lilliput.h"
#include "quarch_benchmark.h"

extern std::mt19937_64 RNG;

void test_circuit_build();
void test_detector_build();
void test_qec_build();
void test_decoding_build();
void test_lilliput_build();

#endif
