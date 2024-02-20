/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 *
 *  An implementation of Kubica and Delfosse's
 *  Restriction Decoder for color codes.
 * */

#ifndef DECODER_RESTRICTION_h
#define DECODER_RESTRICTION_h

#include "qontra/decoder/matching_base.h"

namespace qontra {

typedef std::pair<sptr<graph::decoding::vertex_t>, sptr<graph::decoding::vertex_t>> vpair_t;

class RestrictionDecoder : public MatchingBase {
public:
    typedef std::tuple<uint64_t, uint64_t, int, int> assign_t;
    struct component_t {
        std::vector<assign_t> assignments;
        int color;
    };

    RestrictionDecoder(const DetailedStimCircuit& circuit)
        :MatchingBase(circuit, 3)
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
protected:
    std::vector<component_t>
        compute_connected_components(const std::vector<assign_t>&);

    std::set<vpair_t> insert_error_chain_into(
                        std::map<vpair_t, size_t>& incidence_map,
                        sptr<graph::decoding::vertex_t>,
                        sptr<graph::decoding::vertex_t>,
                        int component_color,
                        int matching_color_1,
                        int matching_color_2,
                        bool force_unflagged=false);

    void insert_incident_vertices(
            std::set<sptr<graph::decoding::vertex_t>>&, const std::set<vpair_t>&, int);
    void insert_incident_vertices(
            std::set<sptr<graph::decoding::vertex_t>>&, const std::map<vpair_t, size_t>&, int);

    std::set<sptr<graph::decoding::hyperedge_t>> get_faces(sptr<graph::decoding::vertex_t>);
};


RestrictionDecoder::assign_t cast_assign(Decoder::assign_t, int, int);
vpair_t make_vpair(sptr<graph::decoding::vertex_t>, sptr<graph::decoding::vertex_t>);

void intersect_with_boundary(
            std::set<vpair_t>&,
            stim::simd_bits_range_ref<SIMD_WIDTH>,
            fp_t& probability,
            sptr<graph::decoding::hyperedge_t>,
            sptr<graph::decoding::vertex_t> incident_vertex);


}   // qontra

#include "restriction.inl"

#endif  // DECODER_RESTRICTION_h
