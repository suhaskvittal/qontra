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
    help += "\tCode name (--code), do not have to use --input and --output if used\n";
    help += "\t\tAssumes that there is a ../graphs/tanner/*.txt file and ../data folder";
    help += "\tInput tanner graph file (--input)\n";
    help += "\tOutput folder (--output)\n";
    help += "\tVerbose (-v)\n";

    std::string code_name;
    std::string file_in;
    std::string folder_out;

    if (parser.get_string("code", code_name)) {
        file_in = std::string("../graphs/tanner/") + code_name + std::string(".txt");
        folder_out = std::string("../data/") + code_name;
    } else {
        parser.get_string("input", file_in);
        parser.get_string("output", folder_out);
    }

    bool verbose = parser.option_set("v");

    // Read the input file.
    std::ifstream fin(file_in);
    io::callback_t<TannerGraph> cb = [] (TannerGraph& g, std::string ln) {
        io::update_tanner_graph(g, ln);
    };
    TannerGraph* graph = new TannerGraph(create_graph_from_file(fin, cb));

    // Print out adjacency list for each parity qubit.
    std::cout << "Z checks:\n";
    for (auto v : graph->get_vertices_by_type(tanner::vertex_t::ZPARITY)) {
        std::cout << "\tZ" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph->get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }
    std::cout << "X checks:\n";
    for (auto v : graph->get_vertices_by_type(tanner::vertex_t::XPARITY)) {
        std::cout << "\tX" << (v->id & 0x3fff'ffff) << " =";
        for (auto w : graph->get_neighbors(v)) {
            std::cout << " D" << w->id;
        }
        std::cout << "\n";
    }
    // Define the cost function here.
    compiler::cost_t cf = [] (compiler::ir_t* ir)
    {
        // Minimize overall connectivity.
        fp_t connectivity = ir->arch->get_mean_connectivity();
        fp_t size = 0.05 * ((fp_t)ir->arch->get_vertices().size());
        fp_t op_cnt = 0.005 * ir->schedule.size();

        fp_t op_depth = 0.0;
        if (ir->dependency_graph != nullptr) {
            op_depth = 0.1 * ir->dependency_graph->get_depth();
        }

        return connectivity + size + op_cnt + op_depth;
    };
    // Define any constraints here.
    compiler::constraint_t con;
    con.max_mean_connectivity = 2.9;
    con.max_connectivity = 4;

    // Declare compiler and run it.
    Compiler compiler(con, cf);
    compiler.params.verbose = verbose;

    compiler::ir_t* res = compiler.run(graph);
    // Print out some simple stats and write to the output folder.
    std::cout << "Cost = " << cf(res) << ", valid = " << res->valid << "\n";
    std::cout << "Number of qubits = " << res->arch->get_vertices().size() << "\n";
    std::cout << "Connectivity = " << res->arch->get_mean_connectivity() << "\n";
    std::cout << "Number of ops = " << res->schedule.size() << "\n";
    std::cout << "Schedule depth = " << res->dependency_graph->get_depth() << "\n";

    write_ir_to_folder(res, folder_out);

    delete res;
    delete graph;

    return 0;
}
