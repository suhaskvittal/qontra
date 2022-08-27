/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef STIM_TEST_h
#define STIM_TEST_h

#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>

#include <stdint.h>

#include <stim.h>

#include "defs.h"
#include "decoding_graph.h"
#include "mwpm_decoder.h"
#include "lilliput.h"
#include "benchmark.h"

extern std::mt19937_64 RNG;

void test_circuit_build();
void test_detector_build();
void test_qec_build();
void test_decoding_build();
void test_lilliput_build();

#endif
