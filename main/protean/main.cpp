/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include <graph/io.h>
#include <graph/tanner_graph.h>
#include <parsing/cmd.h>
#include <protean/network.h>
#include <protean/visualization.h>

using namespace qontra;
using namespace graph;
using namespace protean;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string tanner_graph_file;
    std::string render_output_file;
    std::string data_output_folder;

    if (!pp.get_string("tanner", tanner_graph_file)) return 1;
    if (!pp.get_string("out", data_output_folder)) return 1;

    pp.get_string("render", render_output_file);

    TannerGraph tanner_graph = create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);

    // Make network:
    PhysicalNetwork network(tanner_graph);
    network.join_qubits_with_identical_support();
    // Write data to output folder:
    write_network_to_folder(data_output_folder, network);
#ifdef GRAPHVIZ_ENABLED
    if (render_output_file.size() > 0) {
        // Render the network and save it to a file.
        render_config_t config;
        render_network(render_output_file, network, config);
    }
#endif
    return 0;
}
