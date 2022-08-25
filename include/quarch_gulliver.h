/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef QUARCH_GULLIVER_h
#define QUARCH_GULLIVER_h

#include "quarch_mwpm_decoder.h"
#include "quarch_benchmark.h"

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <set>
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
};

class Gulliver : public MWPMDecoder {
public:
    Gulliver(const stim::Circuit, const GulliverParams&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;
    // Statistics on MWPM vs BFU usage.
    uint32_t n_bfu_accesses;
    uint32_t n_mwpm_accesses;
private:
    // BFUResult holds the matching, the cost, and the cycle time.
    typedef std::tuple<std::map<uint,uint>, fp_t, uint32_t> BFUResult;

    // Recursively examine all possible matchings given a syndrome.
    std::vector<BFUResult> brute_force_matchings(
            const std::vector<uint>&, const BFUResult& running_result);

    uint n_bfu;
    uint32_t n_bfu_cycles_per_add;
    uint bfu_hw_threshold;
    fp_t clock_frequency;
};


#endif
