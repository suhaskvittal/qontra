/*
 *  author: Suhas Vittal
 *  date:   1 September 2023
 * */

#include "mldh/block.h"

namespace qontra {

using namespace graph;

namespace mldh {

Decoder::result_t
BlockDecoder::decode_error(const syndrome_t& syndrome) {
    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    // If the HW is below a given threshold, then there is
    // no reason to spend time computing the blocks. Decode
    // the syndrome directly.
    if (detectors.size() <= config.blocking_threshold) {
        // Update statistics.
        update_blk_statistics(1);
        update_hw_statistics(detectors.size());
        return base_decoder->decode_error(syndrome);
    }
    std::vector<block_t> blocks = get_blocks(detectors);
    // Now, we call our base decoder on each block, and merge
    // the correction together.
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();
    stim::simd_bits corr(obs);
    corr.clear();

    fp_t max_dec_time = 0.0;
    for (const auto& blk : blocks) {
        syndrome_t bs(det+obs);
        bs.clear();

        for (uint x : blk)  {
            bs[x] = 1;
        }
        auto blk_res = base_decoder->decode_error(bs);
        corr ^= blk_res.corr;
        max_dec_time = blk_res.exec_time > max_dec_time ? blk_res.exec_time : max_dec_time;

        // Update block statistics.
        const uint hw = blk.size();
        update_hw_statistics(hw);
    }
    // Update statistics.
    const uint64_t sz = blocks.size();
    update_blk_statistics(sz);

    fp_t t = max_dec_time;

    Decoder::result_t res = {
        t,
        corr,
        is_error(corr, syndrome)
    };
    return res;
}

void
BlockDecoder::update_hw_statistics(uint hw) {
    total_hw_in_block += hw;
    total_hw_sqr_in_block += hw*hw;
    max_hw_in_block = hw > max_hw_in_block ? hw : max_hw_in_block;
    total_blk_hw_above_th += (hw > config.blocking_threshold);
}

void
BlockDecoder::update_blk_statistics(uint sz) {
    total_number_of_blocks += sz;
    total_number_sqr_of_blocks += sz*sz;
    max_number_of_blocks = sz > max_number_of_blocks ? sz : max_number_of_blocks;
    min_number_of_blocks = sz < min_number_of_blocks ? sz : min_number_of_blocks;
    total_shots_evaluated++;
}

std::vector<block_t>
BlockDecoder::get_blocks(std::vector<uint> detectors) {
    // We want to compute the block partitions via a
    // union-find like approach. That is, we grow each
    // block until we cannot anymore.
    std::map<uint, block_t> block_map;
    std::map<uint, uint> root_table;
    for (uint x : detectors) {
        block_map[x] = compute_block_from(x, detectors);
        root_table[x] = x;
    }

    for (uint x : detectors) {
        uint xrt = find(x, root_table);
        block_t& xblk = block_map[xrt];
        // Check if we need to merge anything with this block.
        for (uint y : xblk) {
            uint yrt = find(y, root_table);
            block_t& yblk = block_map[yrt];
            if (yrt == xrt)                 continue;
            if (!block_contains(x, yblk))   continue;
            // Merge yrt's block with xrt.
            merge(xrt, yrt, block_map, root_table);
        }
    }

    std::vector<block_t> blocks;
    for (auto& pair : block_map) {
        uint x = pair.first;
        uint xrt = find(x, root_table);
        // Skip detectors where the detector is
        // not a root.
        if (x != xrt)   continue;
        blocks.push_back(pair.second);
    }
    return blocks;
}

block_t
BlockDecoder::compute_block_from(uint d, std::vector<uint> detectors) {
    auto dv = decoding_graph.get_vertex(d);
    // Get all detectors within distance block_dim of d.
    block_t blk{d};

    fp_t min_edge_weight = std::numeric_limits<fp_t>::max();
    std::map<uint, fp_t> weight_table;
    for (uint x : detectors) {
        if (x == d) continue;
        auto dw = decoding_graph.get_vertex(x);
        auto error_data = decoding_graph.get_error_chain_data(dv, dw);
        if (error_data.chain_length <= block_dim
                && !error_data.error_chain_runs_through_boundary)
        {
            blk.push_back(x);
            min_edge_weight = error_data.weight < min_edge_weight
                                ? error_data.weight : min_edge_weight;
            weight_table[x] = error_data.weight;
        }
    }
    if (config.allow_adaptive_blocks) {
        while (blk.size() > config.cutting_threshold) {
            auto it = std::max_element(blk.begin(), blk.end(),
                        [&] (uint x, uint y) {
                            return weight_table[x] < weight_table[y];
                        });
            if (weight_table[*it] > min_edge_weight + 0.25*block_dim) {
                blk.erase(it);
            } else {
                break;
            }
        }
    }
    return blk;
}

bool
block_contains(uint d, const block_t& blk) {
    return std::find(blk.begin(), blk.end(), d) != blk.end();
}

uint
find(uint x, std::map<uint, uint>& root_table) {
    std::vector<uint> dpath;
    while (x != root_table[x]) {
        dpath.push_back(x);
        x = root_table[x];
    }
    // Now x = root. Update all nodes on dpath.
    for (uint y : dpath)    root_table[y] = x;
    return x;
}

void
merge(uint xrt, uint yrt,
        std::map<uint, block_t>& block_map,
        std::map<uint, uint>& root_table) {
    block_t& xblk = block_map[xrt];
    block_t& yblk = block_map[yrt];

    root_table[yrt] = xrt;
    // Merge blocks.
    for (uint z : yblk) {
        if (!block_contains(z, xblk))   xblk.push_back(z);
    }
    yblk = xblk;
}

}   // mldh
}   // qontra
