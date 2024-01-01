/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#ifndef TESTS_GRAPH_PREBUILT_h
#define TESTS_GRAPH_PREBUILT_h

#include <graph/graph.h>

using namespace qontra;
using namespace graph;
using namespace base;

typedef vertex_t                VERTEX;
typedef edge_t                  EDGE;
typedef Graph<VERTEX, EDGE>     GRAPH;

GRAPH make_complete_graph(size_t);
GRAPH make_bipartite_complete_graph(size_t, size_t);

GRAPH       fast_make_graph_with_k_vertices(size_t k);
sptr<EDGE>  fast_make_edge(sptr<VERTEX>, sptr<VERTEX>);

#include "graph_prebuilt.inl"

#endif  // TESTS_GRAPH_PREBUILT_h
