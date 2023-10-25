/*
 *  author: Suhas Vittal
 *  date:   2 June 2023
 * */

#ifndef GRAPH_IO_h
#define GRAPH_IO_h

#include "defs.h"

#include <fstream>
#include <iostream>
#include <string>

namespace qontra {
namespace graph {

namespace io {
    template <class G_t>
    using callback_t = void (*)(G_t&, std::string); // This call back takes in a line
                                                    // (say from a file) and updates
                                                    // the referenced graph.
}   // io

template <class G_t> G_t
create_graph_from_file(std::string file, io::callback_t<G_t> cb) {
    std::ifstream fin(file);

    G_t graph;
    std::string ln;
    while (std::getline(fin, ln)) {
        cb(graph, ln);
    }
    return graph;
}

template <class G_t> G_t
create_graph_from_string(std::string desc, io::callback_t<G_t> cb) {
    G_t graph;
    
    std::string ln;
    uint prev_index = 0;
    uint index;
    while ((index=desc.find("\n", prev_index)) != std::string::npos) {
        ln = desc.substr(prev_index, index);
        cb(graph, ln);
        prev_index = index+1;
    }
    return graph;
}

}   // graph
}   // qontra

#endif  // GRAPH_IO_h
