/*
 *  author: Suhas Vittal
 *  date:   31 August 2023
 * */

#ifndef MLDH_BLOCK_DECODER_h
#define MLDH_BLOCK_DECODER_h

#include "decoder/mwpm.h"
#include "graph/algorithms/search.h"
#include "graph/decoding_graph.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include <stim.h>
#include <PerfectMatching.h>

namespace qontra {
namespace mldh {

// This is a different approach to decoding,
// where we partition the decoding problem
// into blocks of arbitrary size. Thus, the
// base decoder for a Block Decoder must be
// capable of decoding blocks of arbitrary
// size in space and time.

typedef std::vector<uint> block_t;

class BlockDecoder : public Decoder {
public:
    BlockDecoder(const stim::Circuit& target_circuit,
            Decoder* base,
            const uint block_dim)
        :Decoder(target_circuit),
        base_decoder(base),
        block_dim(block_dim)
    {}

    Decoder::result_t decode_error(const syndrome_t&) override;
    std::string name(void) override { return "Block"; }

    void reset_stats(void) {
        total_ratio_mbhw_thw = 0;
        total_ratio_mbhw_thw_sqr = 0;
        total_ratio_mbhw_thw_max = 0;
        total_shots_above_blk_th = 0;

        total_hw_in_block = 0;
        total_hw_sqr_in_block = 0;
        max_hw_in_block = 0;
        total_blocks = 0;

        total_local_uses = 0;
        total_global_uses = 0;
        total_local_uses_sqr = 0;
        total_global_uses_sqr = 0;
        max_local_uses = 0;
        max_global_uses = 0;
    }


    // Statistics:
    fp_t total_ratio_mbhw_thw=0.0;  // mbhw_thw --> max block hamming weight to total hamming weight.
    fp_t total_ratio_mbhw_thw_sqr=0.0;
    fp_t total_ratio_mbhw_thw_max=0.0;
    uint64_t    total_shots_above_blk_th=0;

    uint64_t    total_hw_in_block=0;
    uint64_t    total_hw_sqr_in_block=0;
    uint64_t    max_hw_in_block=0;
    uint64_t    total_blocks=0;

    uint64_t    total_local_uses=0;
    uint64_t    total_global_uses=0;
    uint64_t    total_local_uses_sqr=0;
    uint64_t    total_global_uses_sqr=0;
    uint64_t    max_local_uses=0;
    uint64_t    max_global_uses=0;


    // Configuration:
    struct {
        bool    auto_to_global_decoder = false;

        uint    blocking_threshold = 14;
        uint    cutting_threshold = 24;
        bool    allow_adaptive_blocks = false;

        fp_t    global_decoder_io_delay = 500;

        struct {
            bool            record_data = false;
            std::ofstream   fout;
        } timing_io;

        struct {
            bool            record_data = false;
            std::ofstream   fout;

            uint64_t        hw_trigger; // If HW > trigger, write data.
            MWPMDecoder*    base;
        } syndrome_io;
    } config;
private:
    typedef std::tuple<int, int, fp_t> blossom_edge_t;

    void update_use_statistics(uint, uint);
    void update_hw_statistics(uint hw);
    void update_blk_ratio_statistics(uint max_blk_hw, uint total_hw);

    void write_timing_data(fp_t);
    void write_syndrome_data(const syndrome_t&, const std::vector<block_t>&, bool is_error);

    bool blossom_subroutine(const std::vector<uint>&, const std::vector<blossom_edge_t>&, Decoder::result_t&);

    std::vector<block_t>    get_blocks(std::vector<uint> detectors, std::vector<blossom_edge_t>&);
    block_t                 compute_block_from(uint d, std::vector<uint> detectors);

    Decoder* base_decoder;
    const uint block_dim;
};

bool block_contains(uint d, const block_t&);

uint
find(uint x, std::map<uint, uint>& root_table);
void
merge(uint xrt, uint yrt,
        std::map<uint, block_t>& block_map,
        std::map<uint, uint>& root_table);


}   // mldh
}   // qontra

#endif  // MLDH_BLOCK_DECODER_h
