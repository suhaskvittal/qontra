/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include <qontra/protean/io.h>
#include <qontra/protean/network.h>
#include <qontra/protean/visualization.h>

#include <qontra/graph/io.h>
#include <qontra/graph/tanner_graph.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    vtils::CmdParser pp(argc, argv);

    std::string tanner_graph_file;
    std::string render_output_folder;
    std::string data_output_folder;

    std::string pass_string = "Jid.Ral.Fla.Ral.Jpa.Ral.(Prx.Con.Jpa.Ral)+.Rlb.";
    std::string layout_engine = "neato";

    uint64_t schedule_rounds = 3;
    uint64_t max_connectivity = 4;

    if (!pp.get_string("tanner", tanner_graph_file)) return 1;
    if (!pp.get_string("out", data_output_folder)) return 1;

    pp.get_string("render", render_output_folder);
    pp.get_string("passes", pass_string);
    pp.get_string("layout", layout_engine);
    
    pp.get_uint64("s-rounds", schedule_rounds);
    pp.get_uint64("max-conn", max_connectivity);
    bool schedule_is_mx = pp.option_set("mx");

    bool verbose = pp.option_set("verbose") || pp.option_set("v");

    TannerGraph tanner_graph = create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);
    std::cout << "Data Qubits = " 
            << tanner_graph.get_vertices_by_type(tanner::vertex_t::type::data).size() << "\n"
            << "Checks = " << tanner_graph.get_checks().size() << "\n";

    vtils::safe_create_directory(data_output_folder);
    // Make network:
    PhysicalNetwork network(tanner_graph);
    network.config.max_connectivity = static_cast<size_t>(max_connectivity);
    network.config.rounds = static_cast<size_t>(schedule_rounds);
    network.config.is_memory_x = schedule_is_mx;

    update_network(pass_string, network, verbose);

    std::cout << "Qubits = " << network.n() << "\n"
            << "Couplings = " << network.m() << "\n"
            << "Mean Degree = " << network.get_mean_connectivity() << "\n"
            << "Max Degree = " << network.get_max_connectivity() << "\n"
            << "Thickness = " << network.get_thickness() << "\n";
    // Write data to output folder:
    write_network_to_folder(data_output_folder, network);
#ifdef GRAPHVIZ_ENABLED
    if (render_output_folder.size() > 0) {
        vtils::safe_create_directory(render_output_folder);
        // Render the network and save it to a file.
        render_config_t rconfig;
        rconfig.layout_engine = layout_engine;
        if (layout_engine == "neato") {
            rconfig.do_not_render_out_of_plane_edges = true;
        }
        render_network(render_output_folder + "/all_checks.pdf", network, rconfig);
        
        rconfig.do_not_render_out_of_plane_edges = false;
        render_network_by_check(render_output_folder, "pdf", network, rconfig);
    }
#endif
    return 0;
}
