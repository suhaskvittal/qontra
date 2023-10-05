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
    RestrictionDecoder(stim::Circuit circuit)
        :Decoder(circuit, graph::DecodingGraph::DO_NOT_BUILD),
        c_decoding_graph(graph::to_colored_decoding_graph(circuit))
    {}

    Decoder::result_t   decode_error(const syndrome_t&) override;

    std::string name() override { return "RestrictionDecoder"; }
protected:
    typedef std::tuple<uint64_t, uint64_t, std::string> match_t;
    typedef std::tuple<uint64_t, uint64_t, std::vector<match_t>, std::string> cc_t;

    std::vector<match_t>    blossom_subroutine(const std::vector<uint>&);
    std::vector<cc_t>       compute_connected_components(const std::vector<match_t>&);

    std::set<colored_vertex_t*> get_incident_vertices(std::vector<colored_edge_t*>, std::string);

    ColoredDecodingGraph c_decoding_graph;
};

}   // qontra

#endif  // RESTRICTION_DECODER_h

