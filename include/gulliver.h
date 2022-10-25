/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef GULLIVER_h
#define GULLIVER_h

#include "mwpm_decoder.h"
#include "benchmark.h"
#include "gulliver/simulator.h"

#include <algorithm>
#include <array>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <vector>

namespace qrc {

struct GulliverParams {
    // Fetch width for brute force unit
    uint bfu_fetch_width;             
    // Max Hamming weight to invoke onchip hardware. 
    uint bfu_hw_threshold;
    // Memory Parameters
    uint n_registers;
    // Cache parameters
    uint n_cache_supertags;
    uint n_cache_sets;
    uint n_cache_lines;
    // DRAM parameters
    std::string dram_config_file;
    std::string log_output_directory;

    fp_t main_clock_frequency;   // in Hz
    fp_t dram_clock_frequency;
};

class Gulliver : public MWPMDecoder {
public:
    Gulliver(const stim::Circuit, uint n_detectors_per_round, 
            const GulliverParams&);
    ~Gulliver();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;
    uint64_t dram_cost(void);

    // Statistics
    uint32_t n_total_accesses;
    uint32_t n_logical_failures;
    fp_t max_latency;
    uint64_t max_bfu_cycles;
    uint64_t max_prefetch_cycles;
    uint max_hamming_weight;
    // More statistics are in the simulator.
    gulliver::GulliverSimulator * simulator;
private:
    uint bfu_hw_threshold;
    uint n_rounds;
    fp_t main_clock_frequency;
    fp_t dram_clock_frequency;
    // Delete later.
    dramsim3::MemorySystem * dram;
    std::map<addr_t, bool> * memory_event_table;
};

}  // qrc

#endif
