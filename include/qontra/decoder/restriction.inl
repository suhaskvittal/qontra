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
        std::set<sptr<graph::decoding::vertex_t>>& vertex_set,
        const std::set<vpair_t>& edge_set,
        int color)
{
    for (const vpair_t& e : edge_set) {
        sptr<graph::decoding::vertex_t> v = e.first,
                                         w = e.second;
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline void
RestrictionDecoder::insert_incident_vertices(
        std::set<sptr<graph::decoding::vertex_t>>& vertex_set,
        const std::map<vpair_t, size_t>& edge_map,
        int color)
{
    for (const auto& p : edge_map) {
        const vpair_t& e = p.first;
        sptr<graph::decoding::vertex_t> v = e.first,
                                        w = e.second;
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline RestrictionDecoder::assign_t
cast_assign(Decoder::assign_t x, int c1, int c2) {
    return std::make_tuple(std::get<0>(x), std::get<1>(x), c1, c2);
}

inline vpair_t
make_vpair(sptr<graph::decoding::vertex_t> v, sptr<graph::decoding::vertex_t> w) {
    if (v < w)  return std::make_pair(v, w);
    else        return std::make_pair(w, v);
}

}   // qontra
