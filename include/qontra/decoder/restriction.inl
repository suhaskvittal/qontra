/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

namespace qontra {

inline void
RestrictionDecoder::insert_incident_vertices(
        std::set<decoding::vertex_t>& vertex_set,
        const std::set<vpair_t>& edge_set,
        int color)
{
    for (const vpair_t& e : edge_set) {
        sptr<decoding::vertex_t> v = e.first,
                                 w = e.second;
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline void
RestrictionDecoder::insert_incident_vertices(
        std::set<decoding::vertex_t>& vertex_set,
        const std::map<vpair_t, size_t>& edge_map,
        int color)
{
    for (const auto& p : edge_set) {
        const vpair_t& e = p.first;
        sptr<decoding::vertex_t> v = e.first,
                                 w = e.second;
        if (v->color == color) vertex_set.insert(v);
        if (w->color == color) vertex_set.insert(w);
    }
}

inline RestrictionDecoder::assign_t
cast_assign(Decoder::assign_t x, int c1, int c2) {
    return std::make_tuple(std::get<0>(x), std::get<1>(x), c1, c2);
}

inline void
intersect_with_boundary(
            std::set<vpair_t>& boundary,
            stim::simd_bits_range_ref<SIMD_WIDTH> corr,
            fp_t& probability,
            sptr<decoding::hyperedge_t> e,
            sptr<decoding::vertex_t> v)
{
    // Get edges of the hyperedge.
    for (size_t j = 0; j < e->get_order(); j++) {
        auto x = e->get<decoding::vertex_t>(j);
        for (size_t k = j+1; k < e->get_order(); k++) {
            auto y = e->get<decoding::vertex_t>(k);
            // Make sure that one of x or y is v.
            if (x != v && y != v) continue;
            vpair_t xy = make_vpair(x, y);
            boundary ^= xy;
        }
    }
    for (uint64_t fr : e->frames) corr[fr] ^= 1;
    probability *= e->probability;
}

}   // qontra
