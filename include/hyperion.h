/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef HYPERION_h
#define HYPERION_h

#include "mwpm_decoder.h"
#include "benchmark.h"
#include "hyperion/simulator.h"

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

#include <math.h>

namespace qrc {

struct HyperionParams {
    // Fetch width for brute force unit
    uint bfu_fetch_width;             
    uint bfu_compute_stages;
    uint bfu_priority_queue_size;
    // Memory Parameters
    uint n_registers;
    // uArch Parameters
    bool use_dma;
    bool use_rc;
    bool use_greedy_init;
    // DRAM parameters
    bool use_dram;
    std::string dram_config_file;
    std::string log_output_directory;

    fp_t main_clock_frequency;   // in Hz
    fp_t dram_clock_frequency;
};

class Hyperion : public MWPMDecoder {
public:
    Hyperion(const stim::Circuit, 
            uint n_detectors_per_round, 
            uint32_t weight_filter_cutoff,
            const HyperionParams&);
    ~Hyperion();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;
    uint64_t dram_cost(void);

    // Statistics
    uint64_t n_nonzero_syndromes;
    uint64_t n_hhw_syndromes;
    uint64_t total_bfu_cycles;
    uint64_t total_prefetch_cycles;
    uint64_t total_cycles_to_converge;
    fp_t total_logfilter_savings;
    uint64_t max_bfu_cycles;
    uint64_t max_prefetch_cycles;
    uint64_t max_cycles_to_converge;
    fp_t min_filter_savings;
    uint max_hamming_weight;
    // More statistics are in the simulator.
    hyperion::HyperionSimulator * simulator;
private:
    uint n_rounds;
    fp_t main_clock_frequency;
    fp_t dram_clock_frequency;
    
    MWPMDecoder baseline;
    // Delete later.
    dramsim3::MemorySystem * dram;
    std::map<addr_t, bool> * memory_event_table;
};

}  // qrc

#endif
