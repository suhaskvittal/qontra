/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#ifndef CODEGEN_DRIVER_h
#define CODEGEN_DRIVER_h

#include "codegen/conv.h"
#include "codegen/tiling.h"

namespace cct {

uptr<qg::TannerGraph> get_sample(uint64_t max_qubits, tiling_config_t, int seed=0);
uptr<qg::TannerGraph> get_best_code(uint64_t max_qubits, uint64_t trials, tiling_config_t);

}   // cct

#endif  // CODEGEN_DRIVER_h
