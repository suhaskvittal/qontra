/*
 *  author: Suhas Vittal
 *  date:   31 August 2023
 * */

#ifndef MLDH_BLOCK_DECODER_h
#define MLDH_BLOCK_DECODER_h

#include "decoder/mwpm.h"
#include "graph/algorithms/search.h"
#include "graph/decoding_graph.h"

#include <limits>
#include <map>
#include <vector>

#include <stim.h>

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
        :Decoder(target_circuit, graph::DecodingGraph::Mode::LOW_MEMORY),
        base_decoder(base),
        block_dim(block_dim)
    {}

    Decoder::result_t decode_error(const syndrome_t&) override;

    std::string name(void) override { return "Block"; }

    // Statistics:
    uint64_t    total_number_of_blocks=0;
    uint64_t    total_number_sqr_of_blocks=0;
    uint64_t    max_number_of_blocks=0;
    uint64_t    min_number_of_blocks=std::numeric_limits<uint64_t>::max();
    uint64_t    total_shots_evaluated=0;

    uint64_t    total_hw_in_block=0;
    uint64_t    total_hw_sqr_in_block=0;
    uint64_t    max_hw_in_block=0;

    uint64_t    total_blk_hw_above_10=0;
private:
    std::vector<block_t>    get_blocks(std::vector<uint> detectors);
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
