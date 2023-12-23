/*
 *  author: Suhas Vittal
 *  date:   24 October 2023
 * */

#define PROTEAN_OBS_OPT

#include <graph/graph.h>
#include <graph/io.h>
#include <graph/tanner_graph.h>
#include <instruction.h>
#include <parsing/cmd.h>

#include <protean/scheduler.h>
#include <protean/utils.h>

#include <fstream>
#include <iostream>

using namespace qontra;

void print_tanner_graph(graph::TannerGraph& gr) {
    auto data_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::type::data);
    auto zparity_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::type::zparity);
    auto xparity_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::type::xparity);
    std::cout << "data = " << data_vertices.size() << "\n";
    std::cout << "zparity = " << zparity_vertices.size() << "\n";
    std::cout << "xparity = " << xparity_vertices.size() << "\n";
}

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);
    std::string tanner_graph_file;
    std::string output_file;
    uint rounds;
    bool is_memory_x = false;

    int seed = 0;
    int randomness = 0;
    int max_depth = -1;

    if (!pp.get_string("tanner", tanner_graph_file))    return 1;
    if (!pp.get_uint32("r", rounds)) return 1;
    if (!pp.get_string("out", output_file)) return 1;

    pp.get_int32("max-depth", max_depth);
    pp.get_int32("seed", seed);
    pp.get_int32("rng", randomness);

    is_memory_x = pp.option_set("x");

    // OPTIMIZATION FLAGS
    protean::oracle::EN_OPT_OBS = pp.option_set("oobs");
    pp.get_float("cobs", protean::oracle::C_OPT_OBS);

    graph::TannerGraph tg = graph::create_graph_from_file(tanner_graph_file, &graph::io::update_tanner_graph);
    print_tanner_graph(tg);

    protean::css_code_data_t code_data;
    if (seed < 0) {
        // Do an exhaustive search for the min depth schedule.
        code_data.schedule_depth = std::numeric_limits<uint>::max();
        const int num_checks = tg.get_vertices_by_type(graph::tanner::vertex_t::type::zparity).size()
                                + tg.get_vertices_by_type(graph::tanner::vertex_t::type::xparity).size();

        std::cout << "Searching for best schedule (seed progress):\n";
        std::cout.flush();

        int rr = 1;
        if (!randomness) rr = 5;
        for (int r = 0; r < rr; r++) {
            std::cout << "\t[round " << r << "]\n";
            for (seed = 0; seed < num_checks; seed++) {
                auto tmp = protean::compute_schedule_from_tanner_graph(tg, seed, randomness, max_depth);
                if (tmp.data_qubits.size() == 0) continue;
                if (pp.option_set("smart-flags")) {
                    tmp = protean::make_fault_tolerant_smart(tmp);
                } else {
                    tmp = protean::make_fault_tolerant_simple(tmp);
                }

                const uint d = tmp.schedule_depth, _d = code_data.schedule_depth;
                const uint f = tmp.flag_qubits.size(), _f = code_data.flag_qubits.size();

#define OBJ(d, f)   ((d) + protean::oracle::C_OPT_OBS*(f))
                if (OBJ(d, f) < OBJ(_d, _f)) {
                    std::cout << "\t" << seed << " (depth = " << d << ", flags = " << f << ")\n";
                    std::cout.flush();
                    code_data = tmp;
                    max_depth = d;
                }
            }
        }
    } else {
        code_data = protean::compute_schedule_from_tanner_graph(tg, seed, randomness);
        if (pp.option_set("smart-flags")) {
            code_data = protean::make_fault_tolerant_smart(code_data);
        } else {
            code_data = protean::make_fault_tolerant_simple(code_data);
        }
    }

    code_data.print_schedule(std::cout);
    std::cout << "qubits: " << code_data.data_qubits.size() << " data, "
                    << code_data.parity_qubits.size() << " parity, "
                    << code_data.flag_qubits.size() << " flags.\n";
    schedule_t mxp = write_memory_experiment(code_data, rounds, is_memory_x);

    std::ofstream fout(output_file);
    fout << schedule_to_text(mxp);

    return 0;
}
