/*
 *  author: Suhas Vittal
 *  date:   05 June 2023
 * */

#include "defs.h"

#include "graph/io.h"
#include "protean/compiler.h"
#include "protean/proc3d.h"
#include "protean/tanner_graph.h"

#include <fstream>
#include <iostream>
#include <string>

#define PRINT_V(id)  (id >> 30) << "|" << ((id >> 24) & 0x3f) << "|" << (id & 0x00ff'ffff)

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    const std::string fname(argv[1]);

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
    // Test out the compiler here with a simple cost function.
    compiler::cost_t cf = [] (Processor3D& proc)
    {
        // Minimize overall connectivity.
        fp_t connectivity = proc.get_connectivity();
        fp_t size_score = 0.05 * ((fp_t)proc.get_vertices().size());
        return connectivity + size_score;
    };

    std::vector<compiler::constraint_t> con;
    compiler::constraint_t con_degree = [] (Processor3D& proc)
    {
        uint max_degree = 0;
        for (auto v : proc.get_vertices()) {
            uint deg = proc.get_degree(v);
            if (deg > max_degree)   max_degree = deg;
        }
        return max_degree <= 4;
    };
    con.push_back(con_degree);
    Compiler compiler(con, cf);
    
    Compiler::ir_t res = compiler.run(graph);

    // Check what the final result is like:
    std::cout << "Number of qubits = " << res.arch.get_vertices().size() << "\n";
    std::cout << "Thickness = " << res.arch.get_thickness() << "\n";
    std::cout << "Connectivity = " << res.arch.get_connectivity() << "\n";
    std::cout << "Number of ops = " << res.schedule.size() << "\n";

    std::cout << "Connections:\n";
    for (auto v : res.arch.get_vertices()) {
        std::cout << "\tQubit " << PRINT_V(v->id) << ":\t";
        for (auto w : res.arch.get_neighbors(v)) {
            std::cout << " " << PRINT_V(w->id);
            if (w->is_tsv_junction())   std::cout << "(V)";
        }
        std::cout << "\n";
    }
    std::cout << "Mapping:\n";
    for (auto pair : res.role_to_qubit) {
        std::cout << "\t" << PRINT_V(pair.first->id) << " --> " << PRINT_V(pair.second->id) << "\n";
    }

    return 0;
}
