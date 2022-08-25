/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef QUARCH_GULLIVER_h
#define QUARCH_GULLIVER_h

#include "quarch_mwpm_decoder.h"
#include "quarch_benchmark.h"

#include "memory_system.h"

#include <algorithm>
#include <array>
#include <condition_variable>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <vector>

struct GulliverParams {
    // Number of functional units for brute force unit
    uint n_bfu;             
    // Number of cycles per addition operation in BFU.
    uint32_t n_bfu_cycles_per_add;
    // Max Hamming weight to invoke BFUs 
    uint bfu_hw_threshold;

    fp_t clock_frequency;   // in Hz
    
    // Memory System Parameters (DramSim3)
    std::string dram_config_file;
    std::string log_output_directory;
};

struct BFUResult {
    std::map<uint, uint> matching;
    fp_t matching_weight;
    bool valid;
};

struct MemoryEvent {
    double tick;
    bool valid;
};

class Gulliver : public MWPMDecoder {
public:
    Gulliver(const stim::Circuit, const GulliverParams&);
    ~Gulliver();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;
    // Statistics on MWPM vs BFU usage.
    uint32_t n_total_accesses;
    uint32_t n_mwpm_accesses;

    uint64_t max_dram_required;
protected:
    /* Computes the address given a detector pair. */
    uint64_t to_addr(uint, uint);
    std::pair<uint, uint> from_addr(uint64_t);

    void bound_detector(uint&);
    void unbound_detector(uint&);
private:
    /* Recursively examine all possible matchings given a syndrome. */
    std::vector<BFUResult> 
        brute_force_matchings(const std::vector<uint>&, uint64_t&);

    dramsim3::MemorySystem * main_memory;
    std::map<std::pair<uint, uint>, MemoryEvent> memory_event_table;
    std::mutex lk_mem_event;
    std::condition_variable cv_table_updated;

    uint n_bfu;
    uint32_t n_bfu_cycles_per_add;
    uint bfu_hw_threshold;
    fp_t clock_frequency;
};


#endif
