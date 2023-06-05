/*
 *  author: Suhas Vittal
 *  date:   05 June 2023
 * */

#include "defs.h"

#include "graph/io.h"
#include "protean/tanner_graph.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    const std::string fname = "../graphs/tanner/rotated_surface_code_d3.txt";

    std::ifstream fin(fname);
    io::callback_t<TannerGraph> cb = [] (TannerGraph& g, std::string ln) {
        io::update_tanner_graph(g, ln);
    };
    TannerGraph graph = create_graph_from_file(fin, cb);

    // Print out adjacency list for each parity qubit.
    std::cout << "Z checks:\n";
    for (auto v : graph.get_vertices_by_type(tanner::vertex_t::ZPARITY)) {
        std::cout << "Z" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph.get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }
    std::cout << "X checks:\n";
    for (auto v : graph.get_vertices_by_type(tanner::vertex_t::XPARITY)) {
        std::cout << "X" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph.get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }

    return 0;
}
