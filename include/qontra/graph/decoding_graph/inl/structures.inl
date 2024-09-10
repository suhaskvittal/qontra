/*
 *  author: Suhas Vittal
 *  date:   29 February 2024
 * */

namespace qontra {
namespace graph {
namespace decoding {

inline sptr<vertex_t>
vertex_t::get_base() {
    return base;
}

}   // decoding

template <> inline std::string
print_v<decoding::vertex_t>(sptr<decoding::vertex_t> v) {
    if (v == nullptr) return "_";
    std::string s;
    s += v->is_boundary_vertex ? "B" : std::to_string(v->id);
    if (v->color != COLOR_ANY) {
        s += "(c=" + std::to_string(v->color) + ")";
    }
    return s;
}

}   // graph
}   // qontra
