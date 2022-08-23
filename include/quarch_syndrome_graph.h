/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QUARCH_SYNDROME_GRAPH_h
#define QUARCH_SYNDROME_GRAPH_h

#include <boost/graph/adjacency_list.hpp>

#include <stim.h>

#include "quarch_defs.h"
#include "quarch_decoding_graph.h"

// Need to use Boost properties for MWPM.
// Unfortunate, but oh well.
typedef boost::property<
            boost::edge_weight_t,
            qfp_t,
            boost::no_property> fp_edge_t;
typedef boost::adjacency_list<
                boost::vecS,
                boost::vecS,
                boost::undirected,
                boost::no_property,
                fp_edge_t> mwpm_graph;

class 


#endif
