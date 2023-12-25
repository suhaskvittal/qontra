/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "graph/algorithms/distance.h"

#include <algorithm>

#include <assert.h>

namespace qontra {

inline bool
is_boundary(sptr<decoding::vertex_t> v) {
    return v->id == BOUNDARY_INDEX;
}

inline bool
is_colored_boundary(sptr<decoding::vertex_t> v) {
    return v->id == RED_BOUNDARY_INDEX
            || v->id == BLUE_BOUNDARY_INDEX
            || v->id == GREEN_BOUNDARY_INDEX;
}

inline std::string
int_to_color(int x) {
    return std::string("rgb"[x]);
}

inline int
color_to_int(std::string x) {
    if (x == "r")       return 0;
    else if (x == "g")  return 1;
    else                return 2;
}

inline DecodingGraph::matrix_entry_t
DecodingGraph::get_error_chain_data(sptr<decoding::vertex_t> v1, sptr<decoding::vertex_t> v2) {
    update_state(); 
    return distance_matrix[v1][v2];
}

inline DecodingGraph::matrix_entry_t
get_error_chain_data_from_flagged_graph(sptr<decoding::vertex_t> v1, sptr<decoding::vertex_t> v2) {
    return flagged_decoding_graph->get_error_chain_data(v1, v2);
}

inline poly_t
DecodingGraph::get_error_polynomial() {
    update_state();
    return error_polynomial;
}

inline fp_t
DecodingGraph::get_expected_errors() {
    update_state();
    return expected_errors;
}

inline void
DecodingGraph::build_distance_matrix() {
    distance_matrix = create_distance_matrix(this, 
            // Weight function:
            [] (sptr<decoding::edge_t> e) { return e->edge_weight; }m
            [&] (sptr<decoding::vertex_t> src,
                sptr<decoding::vertex_t> dst,
                const std::map<sptr<decoding::vertex_t>, fp_t>& dist,
                const std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>& prev)
            { 
                return _dijkstra_cb(src, dst);
            });
}

inline face_t
make_face(sptr<decoding::colored_vertex_t> v1,
            sptr<decoding::colored_vertex_t> v2,
            sptr<decoding::colored_vertex_t> v3)
{
    std::array<sptr<decoding::colored_vertex_t>, 3> tmp{ v1, v2, v3 };
    std::sort(tmp.begin(), tmp.end());
    return std::make_tuple(tmp[0], tmp[1], tmp[2]);
}

inline bool
ColoredDecodingGraph::contains_face(const face_t& fc) {
    return face_frame_map.count(fc);
}

inline std::set<uint>
ColoredDecodingGraph::get_face_frame_changes(const face_t& fc) {
    return face_frame_map[fc];
}

inline fp_t
ColoredDecodingGraph::get_face_probability(const face_t& fc) {
    return face_prob_map[fc];
}

inline DecodingGraph&
ColoredDecodingGraph::operator[](std::string cc) {
    assert(restricted_color_map.count(cc));
    return restricted_graphs.at(restricted_color_map.at(cc));
}

inline DecodingGraph&
ColoredDecodingGraph::operator[](const char* cc) {
    return (*this)[std::string(cc)];
}

}   // qontra
