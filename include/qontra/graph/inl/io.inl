/*
 *  author: Suhas Vittal
 *  date:   2 June 2023
 * */

#include <fstream>
#include <iostream>

namespace qontra {
namespace graph {

template <class GRAPH, class FUNC> inline GRAPH
create_graph_from_file(std::string file, FUNC cb) {
    std::ifstream fin(file);

    GRAPH graph;
    std::string ln;
    while (std::getline(fin, ln)) {
        cb(graph, ln);
    }
    return graph;
}

template <class GRAPH, class FUNC> inline GRAPH
create_graph_from_string(std::string desc, FUNC cb) {
    GRAPH graph;
    
    std::string ln;
    size_t prev_index = 0;
    size_t index;
    while ((index=desc.find("\n", prev_index)) != std::string::npos) {
        ln = desc.substr(prev_index, index);
        cb(graph, ln);
        prev_index = index+1;
    }
    return graph;
}

}   // graph
}   // qontra
