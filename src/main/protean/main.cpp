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
    // Memory experiment parameters
    help += "\tRounds for memory experiment (--rounds)\n";
    help += "\tPhysical error rate (--per)\n";
    help += "\tObservable to measure (--obs), format should be obs1_obs2_obs3_...\n";
    help += "\tShots (--shots)\n";
    help += "\tMemory experiment (-x or -z), default is memory Z\n";

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

    uint rounds;
    fp_t per;
    std::string obs_str;
    uint64_t shots;

    if (!parser.get_uint32("rounds", rounds))   goto help_exit;
    if (!parser.get_float("per", per))          goto help_exit;
    if (!parser.get_string("obs", obs_str))     goto help_exit;
    if (!parser.get_uint64("shots", shots))     goto help_exit;

    bool verbose = parser.option_set("v");
    bool is_memory_x = parser.option_set("x");

    // Parse the observable string.
    std::vector<uint> obs;
    size_t pssi = 0, ssi;
    while (true) {
        if ((ssi = obs_str.find("_", pssi)) == std::string::npos) {
            ssi = obs_str.size();
        }
        uint x = std::stoi(obs_str.substr(pssi, ssi));
        obs.push_back(x);
        pssi = ssi+1;
        if (ssi == obs_str.size()) break;
    }

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
        if (ir->dependency_graph == nullptr)    return 0.0;
        // Build error and time tables.
        // Only need to define H, Mrc, R, and CX.
        ErrorTable et;
        TimeTable tt;
        for (auto v : ir->arch->get_vertices()) {
            et.op1q["H"][v->id] = 0.1*per;
            et.op1q["Mrc"][v->id] = 10*per;
            et.op1q["R"][v->id] = 0.1*per;

            tt.op1q["H"][v->id] = 50;
            tt.op1q["Mrc"][v->id] = 600;
            tt.op1q["R"][v->id] = 50;

            tt.t1[v->id] = 0.5 * (1/per) * 1e3;

            uint dg = ir->arch->get_degree(v);
            fp_t ct_e = 0.01 * per * pow(M_E, dg - 2);
            for (auto w : ir->arch->get_neighbors(v)) {
                auto v_w = std::make_pair(v->id, w->id);
                auto w_v = std::make_pair(w->id, v->id);

                et.op2q["CX"][v_w] = per;
                et.op2q["CX"][w_v] = per;

                tt.op2q["CX"][v_w] = 150;
                tt.op2q["CX"][w_v] = 150;

                if (ct_e > et.op2q_crosstalk["CX"][v_w]) {
                    et.op2q_crosstalk["CX"][v_w] = ct_e;
                    et.op2q_crosstalk["CX"][w_v] = ct_e;
                }
            }
        }

        // Run memory experiment using MWPM decoding.
        write_ir_to_folder(ir, "protean/tmp");
        stim::Circuit circuit = build_stim_circuit(
                                    ir,
                                    rounds,
                                    obs,
                                    is_memory_x,
                                    et,
                                    tt);

        decoder::MWPMDecoder mwpm(circuit);

        experiments::G_USE_MPI = false; // Disable when MPI version is ready.
        auto res = memory_experiment(&mwpm, shots);
        return res.logical_error_rate;
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
