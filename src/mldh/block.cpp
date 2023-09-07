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
    const uint orig_hw = detectors.size();
    // If the HW is below a given threshold, then there is
    // no reason to spend time computing the blocks. Decode
    // the syndrome directly.
    if (orig_hw <= config.blocking_threshold) {
        // Update statistics.
        auto res = base_decoder->decode_error(syndrome);
        write_timing_data(res.exec_time);
        return res;
    }
    std::vector<blossom_edge_t> blossom_edges;
    std::vector<block_t> blocks = get_blocks(detectors, blossom_edges);
    // Now, we call our base decoder on each block, and merge
    // the correction together.
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();
    stim::simd_bits corr(obs);
    corr.clear();

    fp_t max_dec_time = 0.0;
    uint max_blk_hw = 0;
    for (const auto& blk : blocks) {
        syndrome_t bs(det+obs);
        bs.clear();

        for (uint x : blk)  {
            bs[x] = 1;
        }
        Decoder::result_t blk_res;
        blk_res = base_decoder->decode_error(bs);
        corr ^= blk_res.corr;
        max_dec_time = blk_res.exec_time > max_dec_time ? blk_res.exec_time : max_dec_time;

        // Update block statistics.
        const uint hw = blk.size();
        update_hw_statistics(hw);

        max_blk_hw = hw > max_blk_hw ? hw : max_blk_hw;
    }
    // Update statistics.
    update_blk_ratio_statistics(max_blk_hw, orig_hw);
    total_number_of_blocks += blocks.size();
    shots_above_th++;

    fp_t t = max_dec_time;

    Decoder::result_t res = {
        t,
        corr,
        is_error(corr, syndrome)
    };

    write_timing_data(t);
    write_syndrome_data(syndrome, blocks, res.is_error);

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
BlockDecoder::update_blk_ratio_statistics(uint max_blk_hw, uint total_hw) {
    fp_t r = ((fp_t) max_blk_hw) / ((fp_t) total_hw);
    
    total_ratio_mbhw_thw += r;
    total_ratio_mbhw_thw_sqr += r*r;
    total_ratio_mbhw_thw_max = r > total_ratio_mbhw_thw_max ? r : total_ratio_mbhw_thw_max;
}

void
BlockDecoder::write_timing_data(fp_t t) {
    if (!config.timing_io.record_data) {
        return;
    }
    // Write 8 bytes for t (convert to uint64_t).
    uint64_t qt = (uint64_t) t;
    config.timing_io.fout.write(reinterpret_cast<char*>(&qt), sizeof(qt));
}

void
BlockDecoder::write_syndrome_data(const syndrome_t& syndrome, const std::vector<block_t>& blks, bool is_error) {
    if (!config.syndrome_io.record_data) {
        return;
    }
    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    if (detectors.size() <= config.syndrome_io.hw_trigger) {
        return;
    }
    // First, map each detector to a block.
    std::map<uint, uint> detector_to_block_num;
    uint k = 0;
    for (auto blk : blks) {
        for (uint d : blk) {
            detector_to_block_num[d] = k;
        }
        k++;
    }
    Decoder::result_t mwpm_res = config.syndrome_io.base->decode_error(syndrome);
    // We will write the matching in high detail.
    std::ofstream& fout = config.syndrome_io.fout;
    fout << "-------------------------------------------\n";
    fout << "Hamming weight: " << detectors.size() << "\n";
    fout << "Blocks: " << blks.size() << ", sizes =";
    for (auto blk : blks) {
        fout << " " << blk.size();
    }
    fout << "\n";
    fout << "Matching:\n";
    for (auto& mate : mwpm_res.error_assignments) {
        uint x = std::get<0>(mate), y = std::get<1>(mate);
        if (x > y)  continue;
        auto vx = decoding_graph.get_vertex(x);
        auto vy = decoding_graph.get_vertex(y);
        auto error_data = decoding_graph.get_error_chain_data(vx, vy);
        
        const uint len = error_data.chain_length;
        const int bx = detector_to_block_num[x];
        const int by = y == graph::BOUNDARY_INDEX ? -1 : detector_to_block_num[y];
        const bool thru_b = error_data.error_chain_runs_through_boundary;
        if (thru_b && by >= 0) {
            fout << "\t" << x << " (b" << bx << ") <-----[B, " << len << "]-----> " 
                << y << " (b" << by << ")\n";
        } else if (by < 0) {
            fout << "\t" << x << " (b" << bx << ") <-----[" << len << "]-----> B\n";
        } else {
            fout << "\t" << x << " (b" << bx << ") <-----[" << len << "]-----> "
                << y << " (b" << by << ")";
            if (bx != by)   fout << " XX";
            fout << "\n";
        }
        // Print out statistics for x's and y's other options.
        std::map<uint, uint> x_options, y_options;
        for (uint z : detectors) {
            if (x == z || y == z) {
                continue;
            }
            auto vz = decoding_graph.get_vertex(z);
            auto xzdata = decoding_graph.get_error_chain_data(vx, vz);
            const uint xzlen = xzdata.chain_length;
            x_options[xzlen]++;
            if (y != graph::BOUNDARY_INDEX) {
                auto yzdata = decoding_graph.get_error_chain_data(vy, vz);
                const uint yzlen = yzdata.chain_length;
                y_options[yzlen]++;
            }
        }
        fout << "\t\tother options for " << x << ":";
        for (auto pair : x_options) {
            fout << " L" << pair.first << " (" << pair.second << ")";
        }
        fout << "\n";
        if (y != graph::BOUNDARY_INDEX) {
            fout << "\t\tother options for " << y << ":";
            for (auto pair : y_options) {
                fout << " L" << pair.first << " (" << pair.second << ")";
            }
            fout << "\n";
        }
    }
    fout << "Is error && MWPM had no error: " << (is_error && !mwpm_res.is_error) << "\n"; 
}

bool
BlockDecoder::blossom_subroutine(const std::vector<uint>& detectors, 
                                const std::vector<blossom_edge_t>& edges,
                                Decoder::result_t& res)
{
    // This is a shortened implementation of the blossom algorithm.
    // We will really only profile the execution part -- we assume that
    // the bookkeeping aspects of the algorithm can be optimized out
    // (i.e. mapping detectors to ids, and computing the correction)
    // as has been done with PyMatching.
    const uint n = detectors.size();
    const uint m = edges.size();
    // As the existing graph is not guaranteed to have a MWPM, we must
    // reduce it to a graph that will have an MWPM.
    // As we are making this investment, we require that m < O(n^2)/4.
    // Map detectors to ids.
    std::map<uint, uint> detector_to_id;
    for (uint i = 0; i < n; i++) {
        detector_to_id[detectors[i]] = i;
    }
    // Add edges.
    PerfectMatching pm(2*n, 2*(m-n)+n);
    pm.options.verbose = false;

    timer.clk_start();
    uint true_m = 0;
    for (auto& e : edges) {
        int x = std::get<0>(e);
        int y = std::get<1>(e);
        fp_t w = std::get<2>(e);
        if (!detector_to_id.count(x) 
            || (!detector_to_id.count(y) && y >= 0)) 
        {
            continue;
        }

        int i = detector_to_id[x];
        int j = y == -1 ? i+n : detector_to_id[y];
        const uint16_t qw = 1000*w;
        pm.AddEdge(i, j, qw);
        if (y >= 0) {
            // Add duplicate edges.
            pm.AddEdge(i+n, j+n, qw);
        }
        true_m += 1 + (y >= 0);
    }
    if (true_m >= (n*(n-1))/4) {
        timer.clk_end();
        return false;
    }
    pm.Solve();
    fp_t t = timer.clk_end();

    const uint n_observables = circuit.count_observables();
    stim::simd_bits corr(n_observables);
    corr.clear();
    for (uint i = 0; i < n; i++) {
        uint j = pm.GetMatch(i);
        if (i > j)  continue;

        uint x = detectors[i];
        uint y = j == i+n ? graph::BOUNDARY_INDEX : detectors[j];
        auto vx = decoding_graph.get_vertex(x);
        auto vy = decoding_graph.get_vertex(y);
        auto error_data = decoding_graph.get_error_chain_data(vx, vy);
        for (uint fr : error_data.frame_changes) {
            corr[fr] ^= 1;
        }
    }
    res.exec_time = t;
    res.corr = corr;
    return true;
}

std::vector<block_t>
BlockDecoder::get_blocks(std::vector<uint> detectors, std::vector<blossom_edge_t>& blossom_edges) {
    // We want to compute the block partitions via a
    // union-find like approach. That is, we grow each
    // block until we cannot anymore.
    std::map<uint, block_t> block_map;
    std::map<uint, uint> root_table;
    for (uint x : detectors) {
        block_map[x] = compute_block_from(x, detectors);
        root_table[x] = x;
        // Only add blossom edges for those within the same block.
        auto vx = decoding_graph.get_vertex(x);
        for (uint y : block_map[x]) {
            if (x >= y) continue;
            auto vy = decoding_graph.get_vertex(y);
            fp_t w = decoding_graph.get_error_chain_data(vx, vy).weight;
            blossom_edge_t e = std::make_tuple((int)x, (int)y, w);
            blossom_edges.push_back(e);
        }
        // Add an additional edge for the boundary.
        auto vb = decoding_graph.get_vertex(graph::BOUNDARY_INDEX);
        auto wb = decoding_graph.get_error_chain_data(vx, vb).weight;
        blossom_edge_t boundary_e = std::make_tuple((int)x, -1, wb);
        blossom_edges.push_back(boundary_e);
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
    // We compute multiple block rings, where ring r
    // corresponds to all detectors within r of d. We
    // return the smallest block ring r such that block
    // ring r+1 (remove r) is empty (aside from d).
    std::vector<block_t> blk_rings(block_dim);
    for (block_t& b : blk_rings)    b.push_back(d);
    std::vector<uint> len_scatter(block_dim, 0);
    uint min_len = std::numeric_limits<uint>::max();

    fp_t min_edge_weight = std::numeric_limits<fp_t>::max();
    for (uint x : detectors) {
        if (x == d) continue;
        auto dw = decoding_graph.get_vertex(x);
        auto error_data = decoding_graph.get_error_chain_data(dv, dw);
        uint len = error_data.chain_length;
        if (len <= block_dim && !error_data.error_chain_runs_through_boundary) {
            // Note len is 1-indexed, whereas best_blk is 0-indexed.
            len--;  // Subtract 1 from len for this reason.
            for (uint i = len; i < block_dim; i++) {
                blk_rings[i].push_back(x);
            }
            len_scatter[len]++;
            min_len = len < min_len ? len : min_len;
        }
    }
    uint best_blk = 0;
    if (min_len < std::numeric_limits<uint>::max()) {
        best_blk = config.allow_adaptive_blocks ? min_len + block_dim/(2*(min_len+1)) : block_dim-1;
    }
    if (best_blk >= block_dim)  best_blk = block_dim-1;
    return blk_rings[best_blk];
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
