/*
 *  author: Suhas Vittal
 *  date:   8 August 2022
 * */

#ifndef LILLIPUT_h
#define LILLIPUT_h

#include "decoder.h"
#include "mwpm_decoder.h"

#include <boost/math/special_functions/binomial.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <utility>
#include <chrono>

// TODO: Integration with USIMM for DRAM simulation.
struct LILLIPUTParams {
    // System Parameters.
    fp_t clock_frequency;  // In Hz
    bool is_dram_based;
    uint64_t max_memory_size;   // In bytes.
    // DRAM Parameters. All times are in cycles.
    // See: https://www.techpowerup.com/articles/overclocking/64
    uint tCAS;      // Delay between READ command
                    // and data being available.
    uint tRC;       // Row cycle time
    uint tRP;       // Precharge time
    uint tRCD;      // Row address to column address delay

    uint n_channels;
    uint n_ranks;
    uint n_banks;
    uint row_size;  // number of bytes in a row
    uint line_size; // number of bytes in a line
    // SRAM Parameters
    uint sram_access_time;  // in cycles.
    // Decoding Parameters.
    uint detectors_per_round;
};

class LILLIPUT : public Decoder {
/*
 *  An LUT decoder. See: https://arxiv.org/abs/2108.06569
 *
 *  P. Das, A. Locharla, C. Jones. LILLIPUT: A Lightweight Low-Latency Lookup-Table Based Decoder for Near-term Quantum Error Correction. MICRO 2021.
 * */
public:
    LILLIPUT(const stim::Circuit&, Decoder*, 
            const LILLIPUTParams&, std::mt19937_64&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    void clear_stats(void) override;
    void clear_lut(void);

    LILLIPUTParams params;

    uint32_t lut_misses;
protected:
    void populate_lut(std::mt19937_64& rng);
    DecoderShotResult get_lut_entry(const std::vector<uint8_t>&);

    std::vector<DecoderShotResult> lut;
    std::vector<uint8_t> is_lut_entry_valid;
    Decoder * backing_decoder_p;
private:
    fp_t lut_memory_overhead;   // in bytes
    uint lut_entry_size;        // in bytes
};

#endif
