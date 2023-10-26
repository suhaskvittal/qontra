/*
 *  author: Suhas Vittal
 *  date:   24 October 2023
 * */

#include <graph/graph.h>
#include <graph/io.h>
#include <graph/tanner_graph.h>
#include <instruction.h>
#include <parsing/cmd.h>
#include <protean/scheduler.h>

using namespace qontra;

void print_tanner_graph(graph::TannerGraph& gr) {
    auto data_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::Type::data);
    auto zparity_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::Type::zparity);
    auto xparity_vertices = gr.get_vertices_by_type(graph::tanner::vertex_t::Type::xparity);
    std::cout << "data = " << data_vertices.size() << "\n";
    std::cout << "zparity = " << zparity_vertices.size() << "\n";
    std::cout << "xparity = " << xparity_vertices.size() << "\n";
}

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);
    std::string tanner_graph_file;
    int start = 0;

    if (!pp.get_string("tanner", tanner_graph_file))    return 1;
    pp.get_int32("seed", start);

    graph::TannerGraph tg = graph::create_graph_from_file(tanner_graph_file, &graph::io::update_tanner_graph);
    print_tanner_graph(tg);
    auto code_data = protean::compute_schedule_from_tanner_graph(tg, start);

    code_data.print_schedule(std::cout);

    schedule_t mxp = write_memory_experiment(code_data, 4, false);
    std::cout << schedule_to_text(mxp) << "\n";

    return 0;
}
