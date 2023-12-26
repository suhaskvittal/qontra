/*
 *  author: Suhas Vittal
 *  date:   22 June 2023
 * */

#ifndef TIMER_h
#define TIMER_h

#include <stdint.h>
#include <time.h>

// This is a utility class for timing software.
// 
// There are two functions:
//  (1) clk_start: records the time it was called.
//  (2) clk_end: records the time elapsed since the last clk_start.

class Timer {
public:
    void clk_start(void);
    uint64_t clk_end(void);
private:
#ifdef __APPLE__
    uint64_t t_start;
#else
    struct timespec t_start;
#endif
};

#include "timer.inl"

#endif  // TIMER_h
