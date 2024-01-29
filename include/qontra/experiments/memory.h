/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 *
 *  Useful memory experiment functions.
 * */

#ifndef QONTRA_EXPERIMENTS_MEMORY_h
#define QONTRA_EXPERIMENTS_MEMORY_h

#include "qontra/decoder.h"
#include "qontra/ext/stim.h"
#include "qontra/tables.h"

#include <string>

namespace qontra {

struct memory_config_t {
    uint64_t    shots           = 0;
    std::string trace_folder    = "DNE";
};

struct memory_result_t {
    fp_t        logical_error_rate;
    fp_t        hw_mean;
    fp_t        hw_std;
    uint64_t    hw_max;
    fp_t        t_mean;
    fp_t        t_std;
    fp_t        t_max;
    // More specific data:
    std::vector<fp_t>   logical_error_rate_by_obs;  // Logical error rate for each observable.
};

DetailedStimCircuit make_circuit(std::string qes_file, fp_t);
DetailedStimCircuit make_circuit(qes::Program<>, fp_t);
memory_result_t     memory_experiment(Decoder*, memory_config_t);

template <class PROLOGUE, class EPILOGUE>
memory_result_t memory_experiment(Decoder*, memory_config_t, PROLOGUE, EPILOGUE);

}   // qontra

#include "memory.inl"

#endif  // QONTRA_EXPERIMENTS_MEMORY_h
