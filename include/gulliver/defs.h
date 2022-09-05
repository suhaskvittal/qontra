/*
 *  author: Suhas Vittal
 *  date:   01 September 2022
 * */

#ifndef GULLIVER_DEFS_h
#define GULLIVER_DEFS_h

#include <stdint.h>

struct GulliverCycles {
    // Data structure for tracking cycles on-chip and
    // cycles in DRAM.
    uint64_t onchip;
    uint64_t dram;

    GulliverCycles()
        :onchip(0), dram(0)
    {}

    GulliverCycles(uint64_t oc, uint64_t dr)
        :onchip(oc), dram(dr)
    {}

    GulliverCycles operator+(const GulliverCycles& other) const {
        return GulliverCycles(onchip + other.onchip, dram + other.dram);
    }

    GulliverCycles& operator+=(const GulliverCycles& other) {
        onchip += other.onchip;
        dram += other.dram;
        return *this;
    }
};

#endif
