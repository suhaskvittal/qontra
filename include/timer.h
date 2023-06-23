/*
 *  author: Suhas Vittal
 *  date:   22 June 2023
 * */

#ifndef TIMER_h
#define TIMER_h

#include "defs.h"

#include <time.h>

namespace qontra {

// This is a utility class for timing software.
// 
// There are two functions:
//  (1) clk_start: records the time it was called.
//  (2) clk_end: records the time elapsed since the last clk_start.

class Timer {
public:
    void clk_start(void) {
#ifdef __APPLE__
        t_start = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
        struct timespec d;
        clock_gettime(CLOCK_MONOTONIC_RAW, &d);
        t_start = d.tv_nsec;
#endif
    }

    uint64_t clk_end(void) {
#ifdef __APPLE__
        t_end = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
        struct timespec d;
        clock_gettime(CLOCK_MONOTONIC_RAW, &d);
        t_end = d.tv_nsec;
#endif
        return t_end - t_start;
    }
private:
#ifdef __APPLE__
    uint64_t    t_start;
#else
    long        t_start;
#endif
};

}   // qontra

#endif  // TIMER_h
