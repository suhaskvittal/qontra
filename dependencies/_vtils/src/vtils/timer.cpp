/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "vtils/timer.h"

namespace vtils {

void 
Timer::clk_start() {
#ifdef __APPLE__
    t_start = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    clock_gettime(CLOCK_MONOTONIC_RAW, &t_start);
#endif
}

uint64_t
Timer::clk_end() {
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

}   // vtils
