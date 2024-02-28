/*
 *  author: Suhas Vittal
 *  date:   28 February 2024
 * */

namespace qontra {
namespace graph {

inline void
EdgeClass::add_vertex(sptr<decoding::vertex_t> v) {
    for (sptr<decoding::hyperedge_t> e : edges) {
        e->endpoints.push_back(std::reinterpret_pointer_cast<void>(v));
    }
}

inline sptr<decoding::hyperedge_t>
EdgeClass::get_representative() const {
    return rep;
}

inline std::vector<sptr<decoding::hyperedge_t>>
EdgeClass::get_edges() const {
    return edges;
}

}   // graph
}   // qontra
