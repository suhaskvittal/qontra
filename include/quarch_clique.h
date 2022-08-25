/*
 *  author: Suhas Vittal
 *  date:   25 August 2022
 * */

#ifndef QUARCH_CLIQUE_h
#define QUARCH_CLIQUE_h

#include "quarch_mwpm_decoder.h"

#include <boost/graph/adjacency_list.hpp>

#include <map>

struct CliqueParams {
    uint n_cycles_AND;  // Number of cycles for an AND gate
    uint n_cycles_NOT;
    uint n_cycles_XOR;

    fp_t clock_frequency;   // in Hz
                        
    uint detectors_per_round;   // For surface code, this is d*d-1.
};

class CliqueDecoder : public MWPMDecoder {
public:
    CliqueDecoder(const stim::Circuit&, const CliqueParams&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    // More statistics.
    uint32_t n_mwpm_accesses;
    uint32_t n_total_accesses;
private:
    std::map<uint, uint> first_round_degree_table;

    uint n_cycles_AND;
    uint n_cycles_NOT;
    uint n_cycles_XOR;

    uint detectors_per_round;

    fp_t clock_frequency;
};

#endif
