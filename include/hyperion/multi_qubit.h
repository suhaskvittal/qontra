/*
 *  author: Suhas Vittal
 *  date:   8 September 2022
 * */

#ifndef HYPERION_MULTI_QUBIT_h
#define HYPERION_MULTI_QUBIT_h

#include "decoding_graph.h"
#include "hyperion.h"
#include "hyperion/simulator.h"

#include <memory_system.h>
#include <stim.h>

#include <deque>
#include <map>
#include <random>
#include <vector>
#include <utility>

#define HMSIM_DEBUG
#define ROWS_PER_QUBIT  4

namespace qrc {
namespace hyperion {

class HyperionMultiQubitSimulator {
public:
    HyperionMultiQubitSimulator(const std::vector<stim::Circuit>&, 
            uint n_decoders, uint n_detectors_per_round,
            const HyperionParams&);
    ~HyperionMultiQubitSimulator();

    void benchmark(uint32_t shots, std::mt19937_64&);
    void reset_stats(void);

    dramsim3::MemorySystem * dram;
    std::vector<HyperionSimulator*> simulators;

    uint32_t n_timeouts;
    uint32_t n_overflows;
    uint32_t n_logical_errors;
    uint32_t max_faults;
    fp_t max_latency;
private:
    void load_simulator(
            HyperionSimulator*, const std::vector<uint8_t>&, uint qubit_id);
    std::vector<uint8_t> get_correction(const std::map<uint, uint>&, uint qubit_id);

    std::map<addr_t, bool> * memory_event_table;
    std::vector<bool> occupied;
    std::vector<stim::Circuit> circuits;
    std::vector<DecodingGraph> graphs;
    std::vector<graph::PathTable<DecodingGraph::Vertex>> path_tables;
    uint n_rounds;
    fp_t main_clock_frequency;
    fp_t dram_clock_frequency;
};

}  // hyperion
}  // qrc

#endif
