/*
 *  author: Suhas Vittal
 *  date:   4 June 2024
 * */

#include <placc/fpn.h>

#include <qontra/graph/io.h>
#include <qontra/graph/tanner_graph.h>

using namespace qontra;
using namespace graph;
using namespace placc;

int main(int argc, char* argv[]) {
    std::string tanner_graph_file(argv[1]);

    TannerGraph tgr =
        create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);
    for (int c = 0; c < 3; c++) {
        fp_t lat1 = 0, lat2 = 0;
        FPN fpn(&tgr, c);
        fpn.place_flags();
        fpn.place_widowed_qubits();
        auto sch = fpn.phase_one_schedule(lat1);
        sch = fpn.phase_two_schedule(lat2);

        std::cout << "[ Rc = " << c << " ] Latency: P1=" << lat1 << "ns\tP2=" << lat2 << "ns"
                << std::endl;
        std::cout << sch << std::endl;
    }
    return 0;
}
