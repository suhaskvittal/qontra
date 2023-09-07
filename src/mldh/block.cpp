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
    std::vector<block_t> blocks = get_blocks(detectors);
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
        auto blk_res = base_decoder->decode_error(bs);
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
    }
    fout << "Is error && MWPM had no error: " << (is_error && !mwpm_res.is_error) << "\n"; 
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
            if (yrt == xrt)                 continue;
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
    if (config.allow_adaptive_blocks && blk.size() > config.cutting_threshold) {
        const std::vector<fp_t> scale_factors{0.5, 0.4, 0.3, 0.2, 0.1};
        for (auto f : scale_factors) {
            while (blk.size() > config.blocking_threshold) {
                auto it = std::max_element(blk.begin(), blk.end(),
                            [&] (uint x, uint y) {
                                return weight_table[x] < weight_table[y];
                            });
                if (*it == d)   goto adaptive_end;
                if (weight_table[*it] > min_edge_weight*block_dim*f) {
                    blk.erase(it);
                } else {
                    break;
                }
            }
        }
    }
adaptive_end:
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
