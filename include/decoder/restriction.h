/* author: Suhas Vittal
 *  date:   19 July 2023
 *
 *  An implementation of Kubica and Delfosse's
 *  Restriction Decoder for color codes.
 * */

#ifndef RESTRICTION_DECODER_h
#define RESTRICTION_DECODER_h

#include "decoder.h"

#include <PerfectMatching.h>

#include <assert.h>

namespace qontra {

class RestrictionDecoder : public Decoder {
public:
    RestrictionDecoder(DetailedStimCircuit circuit, uint64_t detectors_per_round)
        :Decoder(circuit, graph::DecodingGraph::Mode::DO_NOT_BUILD),
        c_decoding_graph(graph::to_colored_decoding_graph(circuit)),
        flags_are_active(false),
        detectors_per_round(detectors_per_round)
    {}

    Decoder::result_t   decode_error(stim::simd_bits_range_ref) override;

    std::string name() override { return "RestrictionDecoder"; }
protected:
    typedef std::tuple<uint64_t, uint64_t, std::string> match_t;
    typedef std::tuple<uint64_t, uint64_t, std::vector<match_t>, std::string> cc_t;

    template <class VPTR> graph::DecodingGraph::matrix_entry_t
    get_error_chain_data(VPTR x, VPTR y, std::string rc) {
        if (flags_are_active)   return c_decoding_graph[rc].get_error_chain_data_from_flagged_graph(x, y);
        else                    return c_decoding_graph[rc].get_error_chain_data(x, y);
    }

    match_t 
    make_match(uint64_t x, uint64_t y, std::string color) {
        if (x > y) std::swap(x, y);
        return std::make_tuple(x, y, color);
    }

    std::vector<match_t>    blossom_subroutine(const std::vector<uint>&);
    std::vector<cc_t>       compute_connected_components(const std::vector<match_t>&);

    void    insert_error_chain_into(std::set<sptr<graph::decoding::colored_edge_t>>&,
                                    std::map<sptr<graph::decoding::colored_edge_t>, uint>& incidence_map,
                                    std::string component_color,
                                    sptr<graph::decoding::colored_vertex_t> src,
                                    sptr<graph::decoding::colored_vertex_t> dst,
                                    std::string r_color);

    sptr<graph::decoding::colored_vertex_t> flatten(sptr<graph::decoding::colored_vertex_t>);

    std::set<graph::face_t> get_incident_faces(sptr<graph::decoding::colored_vertex_t>);
    stim::simd_bits         get_correction_for_face(graph::face_t);

    // Switches out the boundary connected to other. Used in case we can't make faces with the boundary.
    void switch_out_boundaries(
                    std::set<sptr<graph::decoding::colored_edge_t>>&,
                    std::map<sptr<graph::decoding::colored_edge_t>, uint>& incidence_map);
    // Changes boundary edges to non-boundary edges/
    void remap_boundary_edges(
                    std::set<sptr<graph::decoding::colored_edge_t>>&,
                    std::map<sptr<graph::decoding::colored_edge_t>, uint>& incidence_map);

    graph::ColoredDecodingGraph c_decoding_graph;

    bool            flags_are_active;
    const uint64_t  detectors_per_round;
};

}   // qontra

#endif  // RESTRICTION_DECODER_h

