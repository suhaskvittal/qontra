/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#ifndef TESTS_GRAPH_PREBUILT
#define TESTS_GRAPH_PREBUILT

#include <graph/graph.h>

{   // Typedef block:

using namespace qontra;
using namespace graph;
using namespace base;

typedef vertex_t                VERTEX;
typedef edge_t                  EDGE;
typedef Graph<VERTEX, EDGE>     GRAPH;

}   // Typedefs end.

template <size_t N>
GRAPH make_complete_graph(void);

GRAPH       fast_make_graph_with_k_vertices(size_t k);
sptr<EDGE>  fast_make_edge(sptr<VERTEX>, sptr<VERTEX>);

#include "graph_prebuilt.inl"

#endif  // TESTS_GRAPH_PREBUILT
