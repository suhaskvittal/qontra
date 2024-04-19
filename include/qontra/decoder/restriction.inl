/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#include <vtils/set_algebra.h>

namespace qontra {

inline bool
face_t::operator==(const face_t& other) const {
    return vertices == other.vertices;
}

inline bool
face_t::operator<(const face_t& other) const {
    return vertices < other.vertices;
}

inline void
RestrictionDecoder::insert_incident_vertices(
        std::set<sptr<gd::vertex_t>>& vertex_set,
        const std::set<vpair_t>& edge_set,
        int color)
{
    for (const auto& [v, w] : edge_set) {
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline void
RestrictionDecoder::insert_incident_vertices(
        std::set<sptr<gd::vertex_t>>& vertex_set,
        const std::map<vpair_t, size_t>& edge_map,
        int color)
{
    for (const auto& [e, cnt] : edge_map) {
        const auto& [v, w] = e;
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline vpair_t
make_vpair(sptr<gd::vertex_t> v, sptr<gd::vertex_t> w) {
    if (v < w)  return std::make_pair(v, w);
    else        return std::make_pair(w, v);
}

inline int
color_plus_offset(int c, int off) {
    return (c + off) % 3;
}

}   // qontra
