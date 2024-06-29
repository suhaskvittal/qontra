/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#include "qontra/decoder/restriction.h"

namespace qontra {

void
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

void
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


}   // qontra
