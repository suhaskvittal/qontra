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

struct GulliverParams {
    // Fetch width for brute force unit
    uint bfu_fetch_width;             
    // Max Hamming weight to invoke onchip hardware. 
    uint bfu_hw_threshold;
    // Memory Parameters
    uint n_registers;
    // DRAM parameters
    std::string dram_config_file;
    std::string log_output_directory;

    fp_t main_clock_frequency;   // in Hz
};

class Gulliver : public MWPMDecoder {
public:
    Gulliver(const stim::Circuit, const GulliverParams&);
    ~Gulliver();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;
    uint64_t dram_cost(void);

    // Statistics
    uint32_t n_total_accesses;
    uint32_t n_mwpm_accesses;
    fp_t max_bfu_latency;
    uint64_t max_cycles;
    // More statistics are in the simulator.
    GulliverSimulator * simulator;
private:
    uint bfu_hw_threshold;
    fp_t main_clock_frequency;
    // Delete later.
    dramsim3::MemorySystem * dram;
    std::map<std::pair<uint, uint>, bool> * memory_event_table;
};


#endif
