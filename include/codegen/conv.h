/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#ifndef CODEGEN_CONV_h
#define CODEGEN_CONV_h

#include "codegen/tiling.h"

#include <qontra/graph/tanner_graph.h>
#include <vtils/mat2.h>

namespace qg=qontra::graph;

namespace cct {

uptr<qg::TannerGraph> to_tanner_graph(uptr<TilingGraph>&);
uptr<qg::TannerGraph> to_tanner_graph(const vtils::Mat2&);

vtils::Mat2 to_parity_matrix(uptr<qg::TannerGraph>&);

}   // cct

#endif // CODEGEN_CONV_h
