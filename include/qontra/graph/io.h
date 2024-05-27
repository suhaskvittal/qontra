/*
 *  author: Suhas Vittal
 *  date:   2 June 2023
 * */

#ifndef GRAPH_IO_h
#define GRAPH_IO_h

#include "qontra/defs.h"

#include <string>

namespace qontra {
namespace graph {

template <class GRAPH, class FUNC> inline GRAPH
create_graph_from_file(std::string file, FUNC cb);

template <class GRAPH, class FUNC> inline GRAPH
create_graph_from_string(std::string desc, FUNC cb);

}   // graph
}   // qontra

#include "inl/io.inl"

#endif  // GRAPH_IO_h
