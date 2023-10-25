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

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);
    std::string tanner_graph_file;

    if (!pp.get_string("tanner", tanner_graph_file))    return 1;

    graph::TannerGraph tg = graph::create_graph_from_file(tanner_graph_file, &graph::io::update_tanner_graph);
    auto code_data = protean::compute_schedule_from_tanner_graph(tg);

    code_data.print_schedule(std::cout);

    return 0;
}
