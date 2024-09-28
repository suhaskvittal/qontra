/*
 *  author: Suhas Vittal
 *  date:   9 August 2024
 * */

#ifndef CODEGEN_CONVERT_h
#define CODEGEN_CONVERT_h

#include "codegen/tiling.h"

#include <qontra/graph/tanner_graph.h>

#include <vtils/mat2.h>

namespace qg=qontra::graph;

namespace cgen {

typedef std::tuple<
    sptr<check_t>, int, // Check, side
    sptr<check_t>, int,
    sptr<check_t>, int >
    ctri_t;  // Check triangle, corresponds to a data qubit.

inline std::string print_ctri(ctri_t t) {
    auto& [x,sx,y,sy,z,sz] = t;
    std::string s = "< " + print_check(x) + "." + std::to_string(sx)
                    + " " + print_check(y) + "." + std::to_string(sy)
                    + " " + print_check(z) + "." + std::to_string(sz)
                    + " >";
    return s;
}

ctri_t make_ctri(sptr<check_t>, int, sptr<check_t>, int, sptr<check_t>, int);

qg::TannerGraph to_tanner_graph(const Tiling&);

vtils::Mat2 to_matrix(const Tiling&);
qg::TannerGraph to_tanner_graph(const vtils::Mat2&);

}   // cgen

#endif  // CODEGEN_CONVERT_h
