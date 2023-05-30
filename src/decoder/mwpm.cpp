/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "decoder/mwpm.h"

namespace qontra {
namespace decoder {

Decoder::result_t
MWPMDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    clk_start();

    // Count number of detectors.
    std::vector<uint> detectors;
    for (uint i = 0; i < n_detectors; i++) {
        if (syndrome[i]) {
            detectors.push_back(i);
        }
    }

    if (detectors.size() & 0x1) {
        // Add boundary to matching graph.
        detectors.push_back(graph::BOUNDARY_INDEX);
    }
    // Build Blossom V instance.
    uint n_vertices = detectors.size();
    uint n_edges = n_vertices * (n_vertices + 1) / 2;  // Graph is complete.
    PerfectMatching pm(n_vertices, n_edges);
    pm.options.verbose = false;
    // Add edges.
    for (uint i = 0; i < n_vertices; i++) {
        uint di = detectors[i];
        auto vi = decoding_graph.get_vertex(di);
        for (uint j = i + 1; j < n_vertices; j++) {
            uint dj = detectors[j];
            auto vj = decoding_graph.get_vertex(dj);
            auto error_data = decoding_graph.get_error_chain_data(vi, vj);
            wgt_t edge_weight;
            if (error_data.weight > 1000.0) edge_weight = 10000000;
            else                            edge_weight = MWPM_TO_INT(error_data.weight);
            pm.AddEdge(i, j, edge_weight);
        }
    }
    // Solve instance.
    pm.Solve();
    std::vector<uint8_t> corr(n_observables, 0);
    for (uint i = 0; i < n_vertices; i++) {
        uint j = pm.GetMatch(i);
        uint di = detectors[i];
        uint dj = detectors[j];
        if (di > dj) continue;

        auto vi = decoding_graph.get_vertex(di);
        auto vj = decoding_graph.get_vertex(dj);
        auto error_data = decoding_graph.get_error_chain_data(vi, vj);
        for (auto f : error_data.frame_changes) {
            if (f >= 0) corr[f] ^= 1;
        }
    }
    auto time_taken = (fp_t)clk_end();

    return (Decoder::result_t) { time_taken, corr, is_error(corr, syndrome) };
}

}   // decoder
}   // qontra
