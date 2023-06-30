/*
 *  author: Suhas Vittal
 *  date:   05 June 2023
 * */

#include "decoder/decoder.h"
#include "decoder/mwpm.h"
#include "defs.h"
#include "experiments.h"
#include "graph/io.h"
#include "graph/tanner_graph.h"
#include "parsing/cmd.h"
#include "protean/compiler.h"
#include "protean/proc3d.h"
#include "tables.h"

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

    bool help_requested = parser.option_set("h");
    if (help_requested) {
help_exit:
        std::cout << help;
        return 0;
    }

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
    compiler::cost_t cf = [&] (compiler::ir_t* ir)
    {
        return (fp_t)ir->arch->get_mean_connectivity() 
            + 0.01*ir->arch->get_vertices().size();
    };

    // Define any constraints here.
    compiler::constraint_t con;
    con.max_connectivity = 3;
    con.max_thickness = 1;

    // Declare compiler and run it.
    Compiler compiler(con, cf);
    compiler.params.verbose = verbose;

    schedule_t tmp; // REMOVE
    compiler::ir_t* res = compiler.run(graph, tmp);
    // Print out some simple stats and write to the output folder.
    std::cout << "Cost = " << cf(res) << ", valid = " << res->valid << "\n";
    std::cout << "Number of qubits = " << res->arch->get_vertices().size() << "\n";
    std::cout << "Connectivity = " << res->arch->get_mean_connectivity() << "\n";
    std::cout << "Number of ops = " << res->schedule.size() << "\n";
    std::cout << "Thickness = " << res->arch->get_thickness() << "\n";
//  std::cout << "Schedule depth = " << res->dependency_graph->get_depth() << "\n";

    write_ir_to_folder(res, folder_out);

    delete res;
    delete graph;

    return 0;
}
