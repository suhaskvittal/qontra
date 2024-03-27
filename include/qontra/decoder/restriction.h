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

namespace gd = graph::decoding;

typedef std::pair<sptr<gd::vertex_t>, sptr<gd::vertex_t>> vpair_t;
typedef std::tuple<sptr<gd::vertex_t>, sptr<gd::vertex_t>, sptr<gd::vertex_t>> vtriplet_t;
typedef std::tuple<uint64_t, uint64_t, int, int> c_assign_t;

struct component_t {
    std::vector<c_assign_t> assignments;
    int color;
};

struct face_t {
    std::vector<sptr<gd::vertex_t>> vertices;
    std::set<uint64_t> frames;
    fp_t probability;

    bool operator==(const face_t&) const;
    bool operator<(const face_t&) const;
};

class RestrictionDecoder : public MatchingBase {
public:
    RestrictionDecoder(const DetailedStimCircuit& circuit)
        :MatchingBase(circuit, 3)
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
protected:
    std::vector<component_t>
        compute_connected_components(const std::vector<c_assign_t>&);

    std::set<vpair_t> insert_error_chain_into(
                        std::map<vpair_t, size_t>& incidence_map,
                        sptr<gd::vertex_t>,
                        sptr<gd::vertex_t>,
                        int component_color,
                        int matching_color_1,
                        int matching_color_2,
                        bool force_unflagged);

    // Returns log probability of correction.
    fp_t lifting(
            stim::simd_bits_range_ref<SIMD_WIDTH> corr,
            const std::map<sptr<gd::hyperedge_t>, sptr<gd::hyperedge_t>>& best_rep_map,
            size_t tries=0);

    std::map<sptr<gd::hyperedge_t>, sptr<gd::hyperedge_t>>
        flatten_edge_map(const std::map<sptr<gd::hyperedge_t>, sptr<gd::hyperedge_t>>&);

    void insert_incident_vertices(std::set<sptr<gd::vertex_t>>&, const std::set<vpair_t>&, int);
    void insert_incident_vertices(std::set<sptr<gd::vertex_t>>&, const std::map<vpair_t, size_t>&, int);

    std::set<face_t> get_faces(
        sptr<gd::vertex_t>,
        const std::map<sptr<gd::hyperedge_t>, sptr<gd::hyperedge_t>>& best_rep_map);

    // Data structures:
    // in_cc_map and not_cc_map track the edges present in and outside of the connected components.
    // The specific in_cc edges are in component_edge_sets.
    std::map<vpair_t, size_t> in_cc_map;
    std::map<vpair_t, size_t> not_cc_map;
    std::vector<std::pair<std::set<vpair_t>, int>> component_edge_sets;

    typedef std::tuple<sptr<gd::hyperedge_t>, std::vector<sptr<gd::vertex_t>>, std::map<vpair_t, size_t>&>
        flag_entry_t;

    std::vector<flag_entry_t> triggered_flag_edges;
};

c_assign_t cast_assign(Decoder::assign_t, int, int);
vpair_t make_vpair(sptr<gd::vertex_t>, sptr<gd::vertex_t>);
face_t make_face(sptr<gd::hyperedge_t>);

int color_plus_offset(int, int);

void intersect_with_boundary(
            std::set<vpair_t>&,
            stim::simd_bits_range_ref<SIMD_WIDTH>,
            fp_t& probability,
            const face_t&,
            sptr<gd::vertex_t> incident_vertex);


}   // qontra

#include "restriction.inl"

#endif  // DECODER_RESTRICTION_h
