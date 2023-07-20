/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

#include <assert.h>

namespace qontra {
namespace decoder {

using namespace graph;

Decoder::result_t
RestrictionDecoder::decode_error(const syndrome_t& syndrome) {
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    timer.clk_start();
    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    // Get matching results for each restricted lattice.
    auto r_latRG = decode_restricted_lattice(detectors, Color::red, Color::green);
    auto r_latRB = decode_restricted_lattice(detectors, Color::red, Color::blue);
    // Now, we must merge the matching results.
    //
    // We do this by "jamming" correction strings of different colors at their R vertex
    // meeting point.
    std::map<decoding::vertex_t*, std::array<decoding::vertex_t*, 2>> r_match_table;
    std::vector<Decoder::result_t*> res_array{&r_latRG, &r_latRB};
    for (uint i = 0; i < 2; i++) {
        auto res_p = res_array[i];
        for (auto aa : res_p->error_assignments) {
            uint d1 = std::get<0>(aa);
            uint d2 = std::get<1>(aa);

            auto v1 = decoding_graph.get_vertex(d1);
            auto v2 = decoding_graph.get_vertex(d2);

            Color c1 = detector_to_color[v1];
            Color c2 = detector_to_color[v2];
            if ((v1->id == BOUNDARY_INDEX && c2 == Color::red)
                || (v2->id == BOUNDARY_INDEX && c1 != Color::red)
                || (v2->id != BOUNDARY_INDEX && c2 == Color::red))
            {
                std::swap(v1, v2);
            }
            r_match_table[v1][i] = v2;
        }
    }
    // Now, compute the correction.
    stim::simd_bits corr(n_observables);
    corr.clear();
    for (auto pair : r_match_table) {
        auto rv = pair.first;
        auto mates = pair.second;

        decoding::vertex_t* vj1 = nullptr;
        decoding::vertex_t* vj2 = nullptr;
        for (uint i = 0; i < 2; i++) {
            auto ov = mates[i];
            auto path = decoding_graph.get_error_chain_data(rv, ov).error_chain;
            for (uint j = 2; j < path.size(); j++) {
                auto vx = path[j-1];
                auto vy = path[j];
                auto e = decoding_graph.get_edge(vx, vy);
                for (auto fr : e->frames)   corr[fr] ^= 1;
            }
            if (i == 0) vj1 = path[1];
            else        vj2 = path[1];
        }
        // Jam vj1 and vj2.
        auto jam_frames = decoding_graph.get_error_chain_data(vj1, vj2).frame_changes;
        for (auto fr : jam_frames)  corr[fr] ^= 1;
    }
    fp_t t = (fp_t)timer.clk_end();
    return (Decoder::result_t) { t, corr, is_error(corr, syndrome) };
}

Decoder::result_t
RestrictionDecoder::decode_restricted_lattice(
                            const std::vector<uint>& detectors,
                            Color c1,
                            Color c2)
{
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    stim::simd_bits restricted_syndrome(n_detectors+n_observables);
    restricted_syndrome.clear();
    for (uint d : detectors) {
        auto v = decoding_graph.get_vertex(d);
        Color c = detector_to_color[v];
        restricted_syndrome[d] ^= (c == c1 || c == c2);
    }
    return MWPMDecoder::decode_error(restricted_syndrome);
}

}   // decoder
}   // qontra
