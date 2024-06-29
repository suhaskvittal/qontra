/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

namespace qontra {

inline bool
face_t::operator==(const face_t& other) const {
    return vertices == other.vertices;
}

inline bool
face_t::operator<(const face_t& other) const {
    return vertices < other.vertices;
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
