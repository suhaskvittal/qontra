/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "qontra/graph/algorithms/distance.h"

#include <algorithm>

#include <assert.h>

namespace qontra {
namespace graph {

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

template <> inline std::string
print_v<decoding::vertex_t>(sptr<decoding::vertex_t> v) {
    if (is_boundary(v)) return "B";
    else                return std::to_string(v->id);
}

template <> inline std::string
print_v<decoding::colored_vertex_t>(sptr<decoding::colored_vertex_t> v) {
    if (is_colored_boundary(v)) return "B[" + v->color + "]";
    else                        return std::to_string(v->id) + "[" + v->color + "]";
}

inline std::string
int_to_color(int x) {
    if (x == 0)         return "r";
    else if (x == 1)    return "g";
    else                return "b";
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
    using namespace decoding;
    distance_matrix = create_distance_matrix<vertex_t, edge_t, matrix_entry_t>(this, 
            // Weight function:
            [] (sptr<decoding::edge_t> e) { return e->edge_weight; },
            // Dijkstra callback:
            [&] (sptr<decoding::vertex_t> src,
                sptr<decoding::vertex_t> dst,
                const std::map<sptr<decoding::vertex_t>, fp_t>& dist,
                const std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>& prev)
            { 
                return dijkstra_cb(src, dst, dist, prev);
            });
}

}   // graph
}   // qontra
