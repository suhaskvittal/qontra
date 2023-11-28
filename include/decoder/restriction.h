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
        detectors_per_round(detectors_per_round)
    {}

    Decoder::result_t   decode_error(stim::simd_bits_range_ref) override;

    std::string name() override { return "RestrictionDecoder"; }
protected:
    typedef std::tuple<uint64_t, uint64_t, std::string> match_t;
    typedef std::tuple<uint64_t, uint64_t, std::vector<match_t>, std::string> cc_t;

    match_t 
    make_match(uint64_t x, uint64_t y, std::string color) {
        if (x > y) std::swap(x, y);
        return std::make_tuple(x, y, color);
    }

    graph::DecodingGraph::matrix_entry_t
    get_error_chain_data(uint deta, uint detb, const std::vector<uint>& detectors, std::string rcolor) {
        auto va = c_decoding_graph.get_vertex(deta);
        auto vb = c_decoding_graph.get_vertex(detb);
        return get_error_chain_data(va, vb, detectors, rcolor);
    }

    graph::DecodingGraph::matrix_entry_t
    get_error_chain_data(
            graph::decoding::colored_vertex_t* va,
            graph::decoding::colored_vertex_t* vb,
            const std::vector<uint>& detectors,
            std::string rcolor)
    {
        // First check which flags have been flipped.
        std::set<graph::decoding::vertex_t*> active_flags;
        for (auto df : circuit.flag_detection_events) {
            if (std::find(detectors.begin(), detectors.end(), df) != detectors.end()) {
                active_flags.insert((graph::decoding::vertex_t*) c_decoding_graph.get_vertex(df));
            }
        }
        return c_decoding_graph[rcolor].get_error_chain_data_considering_flags(
                                                        (graph::decoding::vertex_t*)va,
                                                        (graph::decoding::vertex_t*)vb,
                                                        active_flags);
    }

    std::vector<match_t>    blossom_subroutine(const std::vector<uint>&);
    std::vector<cc_t>       compute_connected_components(const std::vector<match_t>&);

    graph::decoding::colored_vertex_t*  flatten(graph::decoding::colored_vertex_t*);

    std::set<graph::face_t> get_incident_faces(graph::decoding::colored_vertex_t*);
    stim::simd_bits         get_correction_for_face(graph::face_t);

    std::set<graph::decoding::colored_vertex_t*>
        get_incident_vertices(std::set<graph::decoding::colored_edge_t*>, std::string);

    graph::ColoredDecodingGraph c_decoding_graph;

    const uint64_t detectors_per_round;
};

}   // qontra

#endif  // RESTRICTION_DECODER_h

