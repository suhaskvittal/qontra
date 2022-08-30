/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef GULLIVER_h
#define GULLIVER_h

#include "mwpm_decoder.h"
#include "benchmark.h"
#include "gulliver/memory_hierarchy.h"

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
    // Number of functional units for brute force unit
    uint n_bfu;             
    // Number of cycles per addition operation in BFU.
    uint32_t n_bfu_cycles_per_add;
    // Max Hamming weight to invoke BFUs 
    uint bfu_hw_threshold;
    fp_t clock_frequency;   // in Hz
    // Cache Parameters
    uint cacheC;
    uint cacheS;
    uint cacheB;
    // DRAM parameters
    std::string dram_config_file;
    std::string log_output_directory;
};

struct BFUResult {
    std::map<uint, uint> matching;
    fp_t matching_weight;
    bool valid;
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

    GulliverCache * cache;

    uint32_t n_total_accesses;
    uint32_t n_mwpm_accesses;
    fp_t max_bfu_latency;
private:
    /* Recursively examine all possible matchings given a syndrome. */
    std::vector<BFUResult> 
        brute_force_matchings(const std::vector<uint>&, uint64_t&);

    uint n_bfu;
    uint32_t n_bfu_cycles_per_add;
    uint bfu_hw_threshold;
    fp_t clock_frequency;

    uint64_t _sram_cost;
};


#endif
