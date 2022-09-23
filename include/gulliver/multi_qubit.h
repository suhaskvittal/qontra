/*
 *  author: Suhas Vittal
 *  date:   8 September 2022
 * */

#ifndef GULLIVER_MULTI_QUBIT_h
#define GULLIVER_MULTI_QUBIT_h

#include "decoding_graph.h"
#include "gulliver.h"
#include "gulliver/simulator.h"

#include <memory_system.h>
#include <stim.h>

#include <deque>
#include <map>
#include <random>
#include <vector>
#include <utility>

#define ROWS_PER_QUBIT  4

namespace qrc {
namespace gulliver {

class GulliverMultiQubitSimulator {
public:
    GulliverMultiQubitSimulator(const std::vector<stim::Circuit>&, 
            uint n_decoders, uint n_detectors_per_round,
            const GulliverParams&);
    ~GulliverMultiQubitSimulator();

    void benchmark(uint32_t shots, std::mt19937_64&);
    void reset_stats(void);

    uint32_t n_timeouts;
    uint32_t n_overflows;
    uint32_t n_uncomputable;
    fp_t max_latency;
private:
    bool load_simulator(
            GulliverSimulator*, const stim::simd_bits_range_ref&, uint qubit_id);

    dramsim3::MemorySystem * dram;
    std::vector<QubitCache*> caches;
    std::map<addr_t, bool> * memory_event_table;

    std::vector<GulliverSimulator*> simulators;
    std::vector<bool> occupied;

    std::vector<stim::Circuit> circuits;
    std::vector<PathTable> path_tables;
    uint n_rounds;
    fp_t main_clock_frequency;
};

}  // gulliver
}  // qrc

#endif
