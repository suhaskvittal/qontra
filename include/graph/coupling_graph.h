/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#ifndef COUPLING_GRAPH_h
#define COUPLING_GRAPH_h

#include "defs.h"
#include "graph/graph.h"

namespace qontra {
namespace graph {

namespace coupling {

struct vertex_t : base::vertex_t {
    int32_t logical_owner;  // Number of logical qubit which contains the qubit.
};

typedef base::edge_t    edge_t;

}   // coupling

typedef Graph<coupling::vertex_t, coupling::edge_t> CouplingGraph;

}   // graph
}   // qontra

#endif  // COUPLING_GRAPH_h
