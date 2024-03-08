/*
 *  author: Suhas Vittal
 *  date:   6 March 2024
 * */

#include <qontra/graph/tanner_graph.h>
#include <qontra/graph/io.h>

#include <stim.h>

using namespace qontra;
using namespace graph;
using namespace tanner;

int main(int argc, char* argv[]) {
    std::string tanner_graph_file(argv[1]);
    TannerGraph gr = 
            create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);

    // Form a Stim circuit for analysis.
    for (size_t i = 0; i <= 1; i++) {
        bool is_x = i & 1;
        int d = gr.compute_code_distance(is_x);
        std::cout << "d" << (is_x ? "x" : "z") << " = " << d << std::endl;
    }
}
