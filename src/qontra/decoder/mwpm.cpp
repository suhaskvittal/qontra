/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "qontra/decoder/mwpm.h"

#include <PerfectMatching.h>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace qontra {

Decoder::result_t
MWPMDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_observables = circuit.count_observables();

    timer.clk_start();
    std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);

    if (detectors.size() & 0x1) {
        // Add boundary to matching graph.
        detectors.push_back(graph::BOUNDARY_INDEX);
    }
    // Build Blossom V instance.
    size_t n_vertices = detectors.size();
    size_t n_edges = n_vertices * (n_vertices + 1) / 2;  // Graph is complete.
    PerfectMatching pm(n_vertices, n_edges);
    pm.options.verbose = false;
    // Add edges.
    for (size_t i = 0; i < n_vertices; i++) {
        uint64_t di = detectors[i];

        auto vi = decoding_graph.get_vertex(di);
        for (size_t j = i + 1; j < n_vertices; j++) {
            uint64_t dj = detectors[j];

            auto vj = decoding_graph.get_vertex(dj);
            auto vi_vj = std::make_pair(vi, vj);
            auto error_data = decoding_graph.get_error_chain_data(vi, vj);

            uint32_t edge_weight;
            if (error_data.weight > 1000.0) edge_weight = 10000000;
            else                            edge_weight = MWPM_TO_INT(error_data.weight);
            pm.AddEdge(i, j, edge_weight);
        }
    }
    // Solve instance.
    pm.Solve();
    stim::simd_bits<SIMD_WIDTH> corr(n_observables);
    corr.clear();

    std::vector<Decoder::assign_t>  error_assignments;
    for (size_t i = 0; i < n_vertices; i++) {
        size_t j = pm.GetMatch(i);
        uint64_t di = detectors[i];
        uint64_t dj = detectors[j];
        if (di > dj) continue;

        auto vi = decoding_graph.get_vertex(di);
        auto vj = decoding_graph.get_vertex(dj);
        auto error_data = decoding_graph.get_error_chain_data(vi, vj);

        stim::simd_bits<SIMD_WIDTH> local_assign(n_observables);
        for (auto f : error_data.frame_changes) {
            if (f >= 0) local_assign[f] ^= 1;
        }
        corr ^= local_assign;
        error_assignments.push_back(std::make_tuple(di, dj, local_assign));
    }
    fp_t t = (fp_t)timer.clk_end();
    Decoder::result_t res = {
        t
        corr,
        error_assignments
    };
    return res;
}

}   // qontra
