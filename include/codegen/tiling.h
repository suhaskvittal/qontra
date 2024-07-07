/*
 *  author: Suhas Vittal
 *  date:   29 June 2024
 * */

#ifndef CODEGEN_TILING_h
#define CODEGEN_TILING_h

#include <qontra/graph.h>

namespace qg=qontra::graph;

namespace cct {

struct shape_t : qg::base::vertex_t {
    size_t sides;
    std::vector<uint64_t> qubits;
    // Neighbors is an array. Position of elements matters. The adjacency list
    // in TilingGraph does not have this property.
    std::vector<sptr<shape_t>> neighbors;
    size_t cnt;
    size_t fptr;
    size_t bptr;
    int color;

    inline uint64_t get_qubit(int i) const {
        return qubits.at( i % sides );
    }

    inline void set_qubit(int i, uint64_t x) {
        qubits[i % sides] = x;
    }

    inline sptr<shape_t> get_neighbor(int i) const {
        return neighbors.at( i % sides );
    }

    inline void set_neighbor(int i, sptr<shape_t> s) {
        cnt += neighbors[i%sides] == nullptr && s != nullptr;
        neighbors[i % sides] = s;
    }
};

struct edge_t : qg::base::edge_t {
    bool is_nonlocal =false;
};

typedef qg::Graph<shape_t, edge_t> TilingGraph;

struct tiling_config_t {
    size_t min_sides    =4;
    size_t max_sides    =12;
    
    int max_oop_allowed_radius =2;

    fp_t base_oop_prob  =0.1;
};

uptr<TilingGraph>
    make_random_tiling(uint64_t max_qubits, tiling_config_t, int seed=0);

}   // cct

namespace qontra {
namespace graph {

template <> inline std::string
print_v(sptr<cct::shape_t> v) {
    if (v == nullptr) return "<nil>";
    return std::to_string(v->id);
}

}
}

#endif  // CODEGEN_TILING_h
