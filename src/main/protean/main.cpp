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

    uint rounds = 3;
    uint64_t shots = 100'000;

    // Define the cost function here.
    compiler::cost_t cf = [&] (compiler::ir_t* ir)
    {
        auto tanner_graph = ir->curr_spec;
        auto arch = ir->arch;
        auto dgr = ir->dependency_graph;
        
        // Generate memory experiment ASM.
        schedule_t mexp;
        //
        // PROLOGUE: initialization.
        //
        const uint n_qubits = arch->get_vertices().size();
        Instruction init_reset;
        init_reset.name = "reset";
        for (uint i = 0; i < n_qubits; i++) {
            init_reset.operands.push_back(i);
        }
        mexp.push_back(init_reset);
        //
        // Syndrome extraction
        //
        uint mctr = 0;
        uint ectr = 0;
        std::vector<tanner::vertex_t*> meas_order;
        for (uint r = 0; r < rounds; r++) {
            uint moffset = mctr;
            for (uint d = 1; d <= dgr->get_depth(); d++) {
                Instruction hsimd = {"h", {}};
                Instruction cxsimd = {"cx", {}};
                Instruction msimd = {"mrc", {mctr}};
                Instruction rsimd = {"reset", {}};
                for (auto v : dgr->get_vertices_at_depth(d)) {
                    auto inst = *(v->inst_p);
                    if (inst.name == "mnrc") {
                        // Convert this to an mrc instruction.
                        for (uint i : inst.get_qubit_operands()) {
                            msimd.operands.push_back(i);
                        }
                        mctr += inst.get_qubit_operands().size();
                        // Update the measurement order.
                        if (r == 0) {
                            auto tv = tanner_graph->get_vertex(v->id);
                            meas_order.push_back(tv);
                        }
                    } else if (inst.name == "h") {
                        for (uint i : inst.operands) {
                            hsimd.operands.push_back(i);
                        }
                    } else if (inst.name == "cx") {
                        for (uint i : inst.operands) {
                            cxsimd.operands.push_back(i);
                        }
                    } else if (inst.name == "reset") {
                        for (uint i : inst.operands) {
                            rsimd.operands.push_back(i);
                        }
                    }
                }    
                std::vector<Instruction*> simd_ops{&hsimd, &cxsimd, &msimd, &rsimd};
                for (auto inst_p : simd_ops) {
                    if (inst_p->get_qubit_operands().size()) mexp.push_back(*inst_p);
                }
            }
            // Add detection events.
            for (uint i = 0; i < meas_order.size(); i++) {
                auto tv = meas_order[i];
                if (tv->qubit_type == tanner::vertex_t::XPARITY)    continue;
                Instruction det;
                det.name = "event";
                det.operands.push_back(ectr++);
                det.operands.push_back(moffset + i);
                if (r > 0) {
                    det.operands.push_back(moffset + i - meas_order.size());
                }
                mexp.push_back(det);
            }
        }
        //
        // Epilogue
        //
        Instruction dqmeas;
        Instruction obs;
        dqmeas.name = "mrc";
        dqmeas.operands.push_back(mctr);

        obs.name = "obs";
        obs.operands.push_back(0);

        std::map<tanner::vertex_t*, uint> data_qubit_meas_order;

        auto data_qubits = tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA);
        for (uint i = 0; i < data_qubits.size(); i++) {
            auto dv = data_qubits[i];
            auto pv = arch->get_vertex(dv->id);
            uint x = ir->qubit_labels[pv];
            dqmeas.operands.push_back(x);
            obs.operands.push_back(mctr + i);

            data_qubit_meas_order[dv] = i;
        }
        mexp.push_back(dqmeas);
        // Add detection events for measurement errors.
        for (uint i = 0; i < meas_order.size(); i++) {
            auto tv = meas_order[i];
            if (tv->qubit_type == tanner::vertex_t::XPARITY)    continue;

            Instruction det;
            det.name = "event";
            det.operands.push_back(ectr++);
            det.operands.push_back(mctr + i - meas_order.size());

            for (auto dv : tanner_graph->get_neighbors(tv)) {
                uint j = data_qubit_meas_order[dv];
                det.operands.push_back(mctr + j);
            }
            mexp.push_back(det);
        }
        mexp.push_back(obs);
        mexp.push_back((Instruction){"done", {}});

        // Now that we have the memory experiment schedule, build the error
        // model. Then, use the Control Simulator to generate a canonical circuit
        // for evaluation.
        //
        // Simple version: uniform model.
        tables::ErrorAndTiming et;
        ErrorTable errors;
        TimeTable timing;
        tables::populate(n_qubits, errors, timing, et);
        stim::Circuit circuit = fast_convert_to_stim(mexp, errors, timing);

        // Build a decoder, and then benchmark it with a memory experiment.
        decoder::MWPMDecoder dec(circuit);
        auto mexp_res = memory_experiment(&dec, shots);
        
        return mexp_res.logical_error_rate;
    };
    experiments::G_USE_MPI = false;

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
