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
#include "sim/control_sim.h"
#include "tables.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tCode name (--code), do not have to use --input and --output if used\n";
    help += "\t\tAssumes that there is a ../graphs/tanner/*.txt file,";
    help += " ../sdl/*.txt file, and ../data folder\n";
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
    std::string sdl_file;

    if (parser.get_string("code", code_name)) {
        file_in = std::string("../graphs/tanner/") + code_name + std::string(".txt");
        folder_out = std::string("../data/") + code_name;
        sdl_file = std::string("../sdl/") + code_name + std::string(".sdl");
    }
    bool verbose = parser.option_set("v");

    // Read the input file.
    std::ifstream fin(file_in);
    io::callback_t<TannerGraph> cb = [] (TannerGraph& g, std::string ln) {
        io::update_tanner_graph(g, ln);
    };
    TannerGraph* graph = new TannerGraph(create_graph_from_file(fin, cb));

    // Define the cost function here.
    compiler::cost_t cf = [&] (compiler::ir_t* ir)
    {
        // Generate memory experiment ASM.
        schedule_t mexp;
        //
        // PROLOGUE: initialization.
        //
        const uint n_qubits = ir->arch->get_vertices().size();
        Instruction init_reset;
        init_reset.name = "reset";
        for (uint i = 0; i < n_qubits; i++) {
            init_reset.operands.push_back(i);
        }
        mexp.push_back(init_reset);
        //
        // Syndrome extraction
        //
        for (auto inst : ir->schedule) {
            if (inst.name == "mnrc") {
                // Convert this to an mrc instruction.
            }
        }
    };

    // Define any constraints here.
    compiler::constraint_t con;
    con.max_connectivity = 4;
    con.max_thickness = 1;

    // Declare compiler and run it.
    Compiler compiler(con, cf);
    compiler.params.verbose = verbose;

    compiler::ir_t* res = compiler.run(graph, sdl_file);
    // Print out some simple stats and write to the output folder.
    std::cout << "Cost = " << cf(res) << ", valid = " << res->valid << "\n";
    std::cout << "Number of qubits = " << res->arch->get_vertices().size() << "\n";
    std::cout << "Number of couplings = " << res->arch->get_edges().size() << "\n";
    std::cout << "Connectivity = " << res->arch->get_mean_connectivity() << "\n";
    std::cout << "Number of ops = " << res->schedule.size() << "\n";
    std::cout << "Thickness = " << res->arch->get_thickness() << "\n";
    std::cout << "Main is planar = " << res->arch->get_main_processor().planar() << "\n";
    std::cout << "Mean coupling length = " << res->arch->get_mean_coupling_length() << "\n";
    std::cout << "Schedule depth = " << res->dependency_graph->get_depth() << "\n";

    write_ir_to_folder(res, folder_out);

    delete res;
    delete graph;

    return 0;
}
