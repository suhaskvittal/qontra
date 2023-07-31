/*
 *  author: Suhas Vittal
 *  date:   05 June 2023
 * */

#include "decoder/decoder.h"
#include "decoder/mwpm.h"
#include "decoder/restriction.h"
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
using namespace decoder;

int main(int argc, char* argv[]) {
    CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tCode name (--code), do not have to use --input and --output if used\n";
    help += "\t\tAssumes that there is a ../graphs/tanner/*.txt file,";
    help += " ../sdl/*.txt file, and ../data folder\n";
    help += "\tMemory experiments rounds (--rounds)\n";
    help += "\tMemory experiments shots (--shots)\n";
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

    uint rounds;
    uint64_t shots;
    if (!parser.get_uint32("rounds", rounds) || !parser.get_uint64("shots", shots)) goto help_exit;

    // Define the cost function here.
    bool is_first_call = true;
    uint run = 0;
    compiler::cost_t cf = [&] (compiler::ir_t* ir)
    {
        auto tanner_graph = ir->curr_spec;
        auto arch = ir->arch;
        auto dgr = ir->dependency_graph;

        const uint n_qubits = arch->get_vertices().size();
        //
        // Define the error model here:
        //
        tables::ErrorAndTiming et;
        /*
        et.e_ro = 0.0;
        et.e_g2q = 0.0;
        et.e_g1q = 0.0;
        et = et * 5;
        */
        ErrorTable errors;
        TimeTable timing;
        tables::populate(n_qubits, errors, timing, et);

        // Generate memory experiment circuit.
        stim::Circuit circuit;
        //
        // PROLOGUE: initialization.
        //
        for (uint i = 0; i < n_qubits; i++) {
            circuit.append_op("R", {i});
            circuit.append_op("X_ERROR", {i}, errors.op1q["reset"][i]);
        }
        //
        // Syndrome extraction
        //
        uint mctr = 0;
        std::vector<std::pair<tanner::vertex_t*, bool>> meas_order;
        // We also want to track colors.
        std::map<uint, uint> measurement_to_color_id;
        std::set<uint> flag_list;
        // Create the main body.
        stim::Circuit body;
        fp_t round_time = 0.0;
        for (uint d = 1; d <= dgr->get_depth(); d++) {
            fp_t layer_time = 0.0;
            for (auto v : dgr->get_vertices_at_depth(d)) {
                auto inst = *(v->inst_p);
                if (inst.name == "mnrc") {
                    auto tv = tanner_graph->get_vertex(inst.metadata.owning_check_id);
                    // Convert this to an mrc instruction.
                    for (uint i = 0; i < inst.operands.size(); i++) {
                        uint j = inst.operands[i];
                        body.append_op("X_ERROR", {j}, errors.op1q["m"][j]);
                        body.append_op("M", {j});
                        // If this is a gauge qubit, then color it differently than
                        // tv.
                        measurement_to_color_id[mctr+i] = 
                            (tv->id % 3 + inst.metadata.is_for_flag) % 3;

                        fp_t t = timing.op1q["m"][j];
                        if (t > layer_time) layer_time = t;
                    }
                    mctr += inst.get_qubit_operands().size();
                    // Update the measurement order.
                    meas_order.push_back(std::make_pair(tv, inst.metadata.is_for_flag)); 
                } else if (inst.name == "h") {
                    for (uint i : inst.operands) {
                        body.append_op("H", {i});
                        body.append_op("DEPOLARIZE1", {i}, errors.op1q["h"][i]);

                        fp_t t = timing.op1q["h"][i];
                        if (t > layer_time) layer_time = t;
                    }
                } else if (inst.name == "cx") {
                    for (uint i = 0; i < inst.operands.size(); i += 2) {
                        uint j1 = inst.operands[i];
                        uint j2 = inst.operands[i+1];
                        auto j1_j2 = std::make_pair(j1, j2);
                        body.append_op("CX", {j1, j2});
                        body.append_op("L_TRANSPORT", {j1, j2}, 
                                        errors.op2q_leakage_transport["cx"][j1_j2]);
                        body.append_op("L_ERROR", {j1, j2}, 
                                        errors.op2q_leakage_injection["cx"][j1_j2]);
                        body.append_op("DEPOLARIZE2", {j1, j2}, 
                                        errors.op2q["cx"][j1_j2]);

                        fp_t t = timing.op2q["cx"][j1_j2];
                        if (t > layer_time) layer_time = t;
                    }
                } else if (inst.name == "reset") {
                    for (uint i : inst.operands) {
                        body.append_op("R", {i});
                        body.append_op("X_ERROR", {i}, errors.op1q["reset"][i]);

                        fp_t t = timing.op1q["reset"][i];
                        if (t > layer_time) layer_time = t;
                    }
                }
            }    
            round_time += layer_time;
        }
        uint ectr = 0;
        for (uint r = 0; r < rounds; r++) {
            // Add detection events.
            for (auto tv : tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA)) {
                auto pv = arch->get_vertex(tv->id);
                uint i = ir->qubit_labels[pv];

                fp_t e1 = (1 - exp(-round_time/timing.t1[i])) * 0.25;
                fp_t e2 = (1 - exp(-round_time/timing.t2[i])) * 0.5 
                            - (1 - exp(-round_time/timing.t1[i])) * 0.25;
                circuit.append_op("X_ERROR", {i}, e1);
                circuit.append_op("Y_ERROR", {i}, e1);
                circuit.append_op("Z_ERROR", {i}, e2);
            }
            circuit += body;
            for (uint i = 0; i < meas_order.size(); i++) {
                auto tv = meas_order[i].first;
                bool is_for_flag = meas_order[i].second;
                if ((tv->qubit_type == tanner::vertex_t::XPARITY) ^ is_for_flag) continue;
                std::cout << "(r = " << r << ") M"
                        << (tv->qubit_type == tanner::vertex_t::XPARITY ? "X" : "Z")
                        << (tv->id & 255) << "("
                        << is_for_flag << ") ---> E" << ectr << "\n";

                uint det1 = (mctr - i) | stim::TARGET_RECORD_BIT;
                if (r == 0) {
                    circuit.append_op("DETECTOR", {det1}, measurement_to_color_id[i]);
                } else {
                    uint det2 = (mctr - i + meas_order.size()) | stim::TARGET_RECORD_BIT;
                    circuit.append_op("DETECTOR", {det1, det2}, measurement_to_color_id[i]);
                }
                if (is_for_flag)    flag_list.insert(ectr);
                ectr++;
            }
        }
        //
        // Epilogue
        //
        // Measure all the data qubits.
        std::map<tanner::vertex_t*, uint> data_qubit_meas_order;
        auto data_qubits = tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA);
        std::vector<uint> epilogue_obs_operands;
        for (uint i = 0; i < data_qubits.size(); i++) {
            auto dv = data_qubits[i];
            auto pv = arch->get_vertex(dv->id);
            uint x = ir->qubit_labels[pv];

            circuit.append_op("X_ERROR", {x}, errors.op1q["m"][x]);
            circuit.append_op("M", {x});

            data_qubit_meas_order[dv] = i;
            epilogue_obs_operands.push_back((data_qubits.size() - i) | stim::TARGET_RECORD_BIT);
        }
        // Add detection events for measurement errors.
        for (uint i = 0; i < meas_order.size(); i++) {
            auto tv = meas_order[i].first;
            bool is_for_flag = meas_order[i].second;
            if (tv->qubit_type == tanner::vertex_t::XPARITY || is_for_flag) {
                continue;
            }

            uint base_det = (mctr - i + data_qubits.size()) | stim::TARGET_RECORD_BIT;
            std::vector<uint> detectors{base_det};
            for (auto dv : tanner_graph->get_neighbors(tv)) {
                uint j = data_qubit_meas_order[dv];
                detectors.push_back((data_qubits.size() - j) | stim::TARGET_RECORD_BIT);
            }
            circuit.append_op("DETECTOR", detectors, measurement_to_color_id[i]);
        }
        circuit.append_op("OBSERVABLE_INCLUDE", epilogue_obs_operands, 0);
        /*
        for (uint i = 0; i < epilogue_obs_operands.size(); i++) {
            circuit.append_op("OBSERVABLE_INCLUDE", {epilogue_obs_operands[i]}, i+1);
        }
        */
        std::ofstream error_model_out(folder_out 
                                        + "/error_model_r" + std::to_string(run) +".stim");
        error_model_out << circuit << "\n";

        // Build a decoder, and then benchmark it with a memory experiment.
        RestrictionDecoder dec(circuit, flag_list);
        experiments::memory_params_t params;
        params.shots = shots;
        auto mexp_res = memory_experiment(&dec, params);

        std::cout << "\tLogical error rate: " << mexp_res.logical_error_rate << "\n";
        
        run++;

        return ir->arch->get_mean_connectivity();
    };
    experiments::G_USE_MPI = false;

    // Define any constraints here.
    compiler::constraint_t con;

    // Declare compiler and run it.
    Compiler compiler(con, cf);
    compiler.params.verbose = verbose;

    compiler::ir_t* res = compiler.run(graph, sdl_file);
    // Print out some simple stats and write to the output folder.
    std::cout << "Cost = " << res->score << ", valid = " << res->valid << "\n";
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
