/*
 *  author: Suhas Vittal
 *  date:   8 September 2022
 * */

#ifndef GULLIVER_MULTI_QUBIT_h
#define GULLIVER_MULTI_QUBIT_h

#include "gulliver/simulator.h"
#include "gulliver.h"

#include <memory_system.h>
#include <stim.h>

#include <vector>

class GulliverMultiQubitSimulator {
public:
    GulliverMultiQubitSimulator(const std::vector<stim::Circuit>&, 
            uint n_decoders, const GulliverParams&);
    ~GulliverMultiQubitSimulator();

    void benchmark(uint32_t shots);

    uint32_t n_logical_failures;
private:
    dramsim3::MemorySystem * dram;
    std::map<addr_t, bool> * memory_event_table;

    std::vector<GulliverSimulator*> simulators;
    std::vector<stim::Circuit> circuits;

    fp_t main_clock_frequency;
};

GulliverMultiQubitSimulator* build_mqsim(const stim::CircuitGenParameters&);

#endif
