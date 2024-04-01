/*
 *  author: Suhas Vittal
 *  date:   31 March 2024
 * */

namespace qontra {
namespace graph {

inline uint64_t
unified_lattice_id(uint64_t id, int c1, int c2) {
    if (c1 > c2) std::swap(c1, c2);
    return id | (static_cast<uint64_t>(c1+1) << 52) | (static_cast<uint64_t>(c2+1) << 56);
}

inline uint64_t
get_unified_lattice_bits_from_id(uint64_t uid, int& c1, int& c2) {
    c1 = static_cast<int>( (uid >> 52) & ((1L << 4)-1) ) - 1;
    c2 = static_cast<int>( uid >> 56 ) - 1;
    return uid >> 52;
}

inline std::string
print_ufl(uint64_t uid) {
    int c1, c2;
    get_unified_lattice_bits_from_id(uid, c1, c2);
    const uint64_t id = uid & ((1L << 52)-1);
    return std::to_string(id) + "(L=" + std::to_string(c1) + "," + std::to_string(c2) + ")";
}

inline std::string
print_ufl(sptr<decoding::vertex_t> v) {
    return print_ufl(v->id);
}

}   // graph
}   // qontra
