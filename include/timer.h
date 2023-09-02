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
        clock_gettime(CLOCK_MONOTONIC_RAW, &t_start);
#endif
    }

    uint64_t clk_end(void) {
#ifdef __APPLE__
        auto t_end = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
        return t_end - t_start;
#else
        struct timespec t_end;
        clock_gettime(CLOCK_MONOTONIC_RAW, &t_end);
        // Compute difference in times.
        const uint64_t B = 1'000'000'000;
        uint64_t time = (B*t_end.tv_sec + t_end.tv_nsec) - (B*t_start.tv_sec + t_start.tv_nsec);
        return time;
#endif
    }
private:
#ifdef __APPLE__
    uint64_t    t_start;
#else
    struct timespec t_start;
#endif
};

}   // qontra

#endif  // TIMER_h
