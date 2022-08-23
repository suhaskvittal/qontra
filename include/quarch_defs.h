/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QUARCH_DEFS_h
#define QUARCH_DEFS_h

#include <limits>
#include <filesystem>

#include <stdint.h>

typedef double      fp_t;   // floating point type
// We define quantized floating point types.
// Used in MWPM.
#if QFP_SIZE==0
typedef int32_t     qfp_t; // quantized floating point type
#elif QFP_SIZE==1
typedef int16_t     qfp_t;
#elif QFP_SIZE==2
typedef int8_t      qfp_t;
#endif

typedef uint16_t uint;
typedef int16_t sint;

#define KB 1024.0
#define MB (KB*1024.0)
#define GB (MB*1024.0)

qfp_t
quantize(fp_t, fp_t fp_max, qfp_t qfp_max);
void 
safe_create_directory(const std::filesystem::path&);

#endif
