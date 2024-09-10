/*
 *  author: Suhas Vittal
 *  date:   22 June 2023
 * */

#ifndef VTILS_TIMER_h
#define VTILS_TIMER_h

#include <stdint.h>
#include <time.h>

// This is a utility class for timing software.
// 
// There are two functions:
//  (1) clk_start: records the time it was called.
//  (2) clk_end: records the time elapsed since the last clk_start.

namespace vtils {

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

}

#endif  // VTILS_TIMER_h
