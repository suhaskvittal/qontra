/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 *
 *  Useful memory experiment functions.
 * */

#ifndef QEC_MEMORY_h
#define QEC_MEMORY_h

#include <qontra/ext/stim.h>

#include <string>

qontra::DetailedStimCircuit make_circuit(std::string qes_file, fp_t);

#include "memory.inl"

#endif  // QEC_MEMORY_h
