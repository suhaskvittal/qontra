/*
 *  author: Suhas Vittal
 *  date:   16 July 2023
 * */

#include "graph/lattice_graph.h"

namespace qontra {
namespace graph {

using namespace lattice;

LatticeGraph::LatticeGraph(void)
    :Graph()
{}

LatticeGraph::LatticeGraph(const LatticeGraph& other)
    :Graph(other),
    x_obs_list(other.x_obs_list),
    z_obs_list(other.z_obs_list)
{}

}   // graph
}   // qontra
