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
    std::set<size_t> occupied_sides;

    int color;
};

struct edge_t : qg::base::edge_t {
    bool is_nonlocal =false;
};

typedef qg::Graph<shape_t, edge_t> TilingGraph;

struct tiling_config_t {
    size_t min_sides    =4;
    size_t max_sides    =12;
    
    int max_oop_allowed_radius =2;

    fp_t base_oop_prob  =0.01;
};

uptr<TilingGraph>
    make_random_tiling(uint64_t max_qubits, tiling_config_t, int seed=0);

}   // cct

#endif  // CODEGEN_TILING_h
