/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <unistd.h>

namespace qontra {

inline void
configure_optimal_batch_size() {
#ifndef L1D_CACHE_LINE_SIZE
#if defined(linux)
    uint64_t cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);  // In bytes.
#else
    uint64_t cache_line_size = 64;
#endif
#else
    uint64_t cache_line_size = L1D_CACHE_LINE_SIZE;
#endif
    experiments::G_SHOTS_PER_BATCH = cache_line_size << 3;   // Need to convert to bits.
}

}   // qontra
