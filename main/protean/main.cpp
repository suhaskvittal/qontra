/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include <protean/io.h>
#include <protean/network.h>
#include <protean/visualization.h>

#include <qontra/graph/io.h>
#include <qontra/graph/tanner_graph.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#ifdef PROTEAN_PERF
#include <vtils/timer.h>
#endif

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    vtils::CmdParser pp(argc, argv, 2);

    std::string tanner_graph_file(argv[1]);
    std::string data_output_folder(argv[2]);

    std::string render_output_folder = data_output_folder + "/render";

    std::string pass_string = "Jid.Ral.Fla.Ral.Jpa.Ral.(Prx.Con.Jpa.Ral)+.Rlb.";
    std::string layout_engine = "neato";

    uint64_t schedule_rounds = 3;
    uint64_t max_connectivity = 4;

    pp.get("passes", pass_string);
    pp.get("layout", layout_engine);
    
    pp.get("s-rounds", schedule_rounds);
    pp.get("max-conn", max_connectivity);

    bool verbose = pp.option_set("verbose") || pp.option_set("v");

    vtils::safe_create_directory(data_output_folder);

#ifdef PROTEAN_PERF
    vtils::Timer timer;
    fp_t t;

    timer.clk_start();
#endif
    
    // Make Network:
    TannerGraph tanner_graph = create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);
    PhysicalNetwork network(&tanner_graph);

    network.config.max_connectivity = static_cast<size_t>(max_connectivity);
    network.config.rounds = static_cast<size_t>(schedule_rounds);
    network.config.force_unopt_flags = pp.option_set("fno-opt-flags");
    network.config.force_xz_flag_merge = pp.option_set("fflag-jid");
    network.config.enable_flag_reduction = pp.option_set("fflag-reduce");
    network.config.enable_proxy_triangle_search = pp.option_set("fproxy-triangle");

    std::cout << "Data Qubits = " 
            << tanner_graph.get_vertices_by_type(tanner::vertex_t::type::data).size() << "\n"
            << "Checks = " << tanner_graph.get_checks().size() << "\n";
    // Run:
    update_network(pass_string, &network, verbose);
    if (pp.option_set("color-checks") && !pp.option_set("skip-schedule")) {
        network.assign_colors_to_checks();
    }
    if (pp.option_set("skip-schedule")) {
        network.config.skip_scheduling = true;
    }
    network.finalize();
    // Write data to output folder:
    write_network_to_folder(data_output_folder, &network);
    std::cout << "Qubits = " << network.n() << "\n"
            << "Couplings = " << network.m() << "\n"
            << "Mean Degree = " << network.get_mean_degree() << "\n"
            << "Max Degree = " << network.get_max_degree() << "\n"
            << "Thickness = " << network.get_thickness() << "\n"
            << "Round Latency = " << network.get_round_latency() << "\n"
            << "Round CNOTs = " << network.get_round_cnots() << "\n";

#ifdef PROTEAN_PERF
    t = timer.clk_end();
    std::cout << "[ protean ] total time for compilation: " << t*1e-9 << "s" << std::endl;
#endif

#ifdef GRAPHVIZ_ENABLED
    if (!pp.option_set("skip-render") && render_output_folder.size() > 0) {
        vtils::safe_create_directory(render_output_folder);
        // Render the network and save it to a file.
        render_config_t rconfig;
        rconfig.layout_engine = layout_engine;
        if (layout_engine == "neato") {
            rconfig.do_not_render_out_of_plane_edges = true;
        }
        render_network(render_output_folder + "/all_checks.pdf", &network, rconfig);
        
        rconfig.do_not_render_out_of_plane_edges = false;
        render_network_by_check(render_output_folder, "pdf", &network, rconfig);
    }
#endif
    return 0;
}
