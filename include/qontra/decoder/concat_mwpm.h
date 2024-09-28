/*
 *  author: Suhas Vittal
 *  daet:   26 May 2024
 * */

#ifndef DECODER_CONCAT_MWPM_h
#define DECODER_CONCAT_MWPM_h

#include "qontra/decoder/matching_base.h"

#include <vtils/bijective_map.h>

namespace qontra {

namespace gd=graph::decoding;

typedef 
    vtils::BijectiveMap<std::pair<sptr<gd::vertex_t>, sptr<gd::vertex_t>>, sptr<gd::vertex_t>>
    evmap_t;

std::pair<sptr<gd::vertex_t>, sptr<gd::vertex_t>>
    make_ev_pair(sptr<gd::vertex_t>, sptr<gd::vertex_t>);

class ConcatMWPMDecoder : public MatchingBase {
public:
    ConcatMWPMDecoder(const DetailedStimCircuit&);

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
    
    fp_t rgb_compute_matching(
            const std::vector<sptr<gd::vertex_t>>&,
            int color,
            stim::simd_bits_range_ref<SIMD_WIDTH>,
            size_t& correction_weight);
private:
    // The parent MatchingBase class already has a DecodingGraph. We will use this
    // graph for restricted lattice matchings. For matchings on the R/G/B-only lattices,
    // we have three DecodingGraphs for that purpose.
    //
    // Note that the RGB-only lattices have one vertex per edge in the color code tiling.
    // That is, each 3-hyperedge can be split into a vertex-vertex and an edge-vertex.
    std::array<uptr<graph::DecodingGraph>, 3> rgb_only_graphs;
    evmap_t edge_vertex_map;
};

}   // qontra

#include "inl/concat_mwpm.inl"

#endif  // DECODER_CONCAT_MWPM_h
