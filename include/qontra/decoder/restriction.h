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

const int COLOR_RED = 0;
const int COLOR_GREEN = 1;

namespace gd = graph::decoding;

typedef std::pair<sptr<gd::vertex_t>, sptr<gd::vertex_t>> vpair_t;
typedef std::tuple<sptr<gd::vertex_t>, sptr<gd::vertex_t>, sptr<gd::vertex_t>> vtriplet_t;

struct component_t {
    std::vector<assign_t> assignments;
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
    // Set chamberland=true to use Chris Chamberland's original version. This is heavily modified to work
    // with hyperbolic color codes as well.
    RestrictionDecoder(const DetailedStimCircuit& circuit, bool chamberland=false)
        :MatchingBase(circuit, 3, false),
        chamberland(chamberland)
    {}

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>) override;
protected:
    virtual std::vector<assign_t> compute_matchings(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome);

    void split_assignment(
            std::vector<assign_t>&,
            const assign_t&,
            const std::map<sptr<gd::hyperedge_t>, size_t>& flag_ctr_map);

    std::vector<component_t> compute_connected_components(const std::vector<assign_t>&);

    std::set<vpair_t> insert_error_chain_into(
                        std::map<vpair_t, size_t>& incidence_map,
                        const std::vector<sptr<gd::vertex_t>>& path,
                        int component_color,
                        int c1,
                        int c2);

    // Returns log probability of correction.
    fp_t lifting(
            stim::simd_bits_range_ref<SIMD_WIDTH> corr,
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome_delta,
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

    bool chamberland;
};

vpair_t make_vpair(sptr<gd::vertex_t>, sptr<gd::vertex_t>);
face_t make_face(sptr<gd::hyperedge_t>);
int color_plus_offset(int, int);

// Helper functions:
void push_back_assignment(std::vector<assign_t>&, const assign_t&);
void erase_from_incidence_map(std::map<vpair_t, size_t>&, const vpair_t&);

void update_correction(
        // To be updated:
        std::map<vpair_t, size_t>& inc_map,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr,
        fp_t& corr_log_pr,
        // Updated by:
        stim::simd_bits_range_ref<SIMD_WIDTH> local_corr,
        const std::set<vpair_t>& local_boundary,
        fp_t local_log_pr);

bool remove_widowed_edges(std::map<vpair_t, size_t>&);

std::set<vpair_t>
find_face_subset_given_cc_map(
        const std::map<vpair_t, size_t>& x_cc_map,
        const std::set<face_t>& faces,
        stim::simd_bits_range_ref<SIMD_WIDTH>,
        fp_t& probability,
        std::vector<face_t>& applied_faces,
        sptr<gd::vertex_t> incident_vertex);

}   // qontra

#include "inl/restriction.inl"

#endif  // DECODER_RESTRICTION_h
