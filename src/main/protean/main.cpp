/*
 *  author: Suhas Vittal
 *  date:   05 June 2023
 * */

#include "defs.h"
#include "graph/io.h"
#include "parsing/cmd.h"
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
    CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tInput tanner graph file (--input)\n";
    help += "\tOutput folder (--output)\n";
    help += "\tVerbose (-v)\n";

    std::string file_in;
    std::string folder_out;

    parser.get_string("input", file_in);
    parser.get_string("output", folder_out);

    bool verbose = parser.option_set("v");

    // Read the input file.
    std::ifstream fin(file_in);
    io::callback_t<TannerGraph> cb = [] (TannerGraph& g, std::string ln) {
        io::update_tanner_graph(g, ln);
    };
    TannerGraph graph = create_graph_from_file(fin, cb);

    // Print out adjacency list for each parity qubit.
    std::cout << "Z checks:\n";
    for (auto v : graph.get_vertices_by_type(tanner::vertex_t::ZPARITY)) {
        std::cout << "\tZ" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph.get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }
    std::cout << "X checks:\n";
    for (auto v : graph.get_vertices_by_type(tanner::vertex_t::XPARITY)) {
        std::cout << "\tX" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph.get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }
    // Define the cost function here.
    compiler::cost_t cf = [] (Processor3D& proc)
    {
        // Minimize overall connectivity.
        fp_t connectivity = proc.get_connectivity();
        fp_t size_score = 0.05 * ((fp_t)proc.get_vertices().size());
        return connectivity + size_score;
    };
    // Define any constraints here.
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

    // Declare compiler and run it.
    Compiler compiler(con, cf);
    Compiler::ir_t res = compiler.run(graph);
    // Print out some simple stats and write to the output folder.
    std::cout << "Number of qubits = " << res.arch.get_vertices().size() << "\n";
    std::cout << "Connectivity = " << res.arch.get_connectivity() << "\n";
    std::cout << "Number of ops = " << res.schedule.size() << "\n";

    write_ir_to_folder(res, folder_out);

    return 0;
}
