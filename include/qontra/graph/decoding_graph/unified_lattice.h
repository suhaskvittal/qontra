/*
 *  author: Suhas Vittal
 *  date:   31 March 2024
 * */

#ifndef QONTRA_DECODING_GRAPH_UNIFIED_LATTICE_h
#define QONTRA_DECODING_GRAPH_UNIFIED_LATTICE_h

#include "qontra/graph/decoding_graph/structures.h"

#include <string>

namespace qontra {
namespace graph {

uint64_t unified_lattice_id(uint64_t, int, int);
uint64_t get_unified_lattice_bits_from_id(uint64_t, int&, int&);

std::string print_ufl(uint64_t);
std::string print_ufl(sptr<decoding::vertex_t>);

}   // graph
}   // qontra

#include "unified_lattice.inl"
 
#endif  // QONTRA_DECODING_GRAPH_UNIFIED_LATTICE_h
