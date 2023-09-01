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
    // TODO: Add timing results.

    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
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
        total_hw_in_block += hw;
        total_hw_sqr_in_block += hw*hw;
        max_hw_in_block = hw > max_hw_in_block ? hw : max_hw_in_block;
        total_blk_hw_above_10 += (hw > 10);
    }
    // Update statistics.
    const uint64_t sz = blocks.size();
    total_number_of_blocks += sz;
    total_number_sqr_of_blocks += sz*sz;
    max_number_of_blocks = sz > max_number_of_blocks ? sz : max_number_of_blocks;
    min_number_of_blocks = sz < min_number_of_blocks ? sz : min_number_of_blocks;
    total_shots_evaluated++;

    fp_t t = max_dec_time;

    Decoder::result_t res = {
        t,
        corr,
        is_error(corr, syndrome)
    };
    return res;
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
            if (yrt == xrt) continue;
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
    block_t blk;

    std::map<decoding::vertex_t*, uint> depth_table;
    depth_table[dv] = 0;
    // We need to perform a custom BFS to compute the elements
    // in the block.
    std::set<decoding::vertex_t*>   visited;
    std::deque<decoding::vertex_t*> bfs;
    bfs.push_back(dv);
    while (bfs.size()) {
        auto v = bfs.front();
        bfs.pop_front();

        if (depth_table[v] <= block_dim
            && block_contains(v->id, detectors)
            && !block_contains(v->id, blk)) 
        {
            blk.push_back(v->id);
        }

        for (auto w : decoding_graph.get_neighbors(v)) {
            if (w->id == BOUNDARY_INDEX)    continue;
            if (visited.count(w))           continue;
            if (!depth_table.count(w)
                || depth_table[v]+1 < depth_table[w])
            {
                depth_table[w] = depth_table[v] + 1;
                bfs.push_back(w);
            }
        }
        visited.insert(v);
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
