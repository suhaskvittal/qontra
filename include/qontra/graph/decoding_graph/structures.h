/*
 *  author: Suhas Vittal
 *  date:   29 February 2024
 * */

#ifndef QONTRA_DECODING_GRAPH_STRUCTURES_h
#define QONTRA_DECODING_GRAPH_STRUCTURES_h

#include "qontra/hypergraph.h"

#include <limits>
#include <unordered_set>

namespace qontra {
namespace graph {

const int COLOR_ANY = -1;
const uint64_t BOUNDARY_FLAG = std::numeric_limits<uint32_t>::max();

namespace decoding {

struct vertex_t : base::vertex_t {
    int             color = COLOR_ANY;
    sptr<vertex_t>  base = nullptr;
    bool            is_boundary_vertex = false;

    sptr<vertex_t> get_base(void);
};

struct edge_t : base::edge_t {
    fp_t probability;
};

struct hyperedge_t : base::hyperedge_t {
    fp_t                probability;
    fp_t                power = 1.0;
    std::unordered_set<uint64_t>  flags;
    std::unordered_set<uint64_t>  frames;
};

}   // decoding
}   // graph
}   // qontra

#include "inl/structures.inl"

#endif  // QONTRA_DECODING_GRAPH_STRUCTURES_h
