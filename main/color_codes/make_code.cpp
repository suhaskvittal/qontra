/*
 *  author: Suhas Vittal
 *  date:   3 July 2024
 * */

#include <codegen/monte_carlo.h>
#include <codegen/convert.h>

#include <vtils/cmd_parse.h>

#include <string>

#include <mpi.h>
#include <stdlib.h>

using namespace qontra;
using namespace graph;
using namespace cgen;
using namespace vtils;

int main(int argc, char* argv[]) {
    int r = atoi(argv[1]);
    int c = atoi(argv[2]);
    int s = atoi(argv[3]);

    CmdParser pp(argc, argv, 3);

    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    s = (s/size) + (rank==0)*(s%size);

    MonteCarloManager m(rank);
    m.config.r = r;
    m.config.c = c;
    pp.get("p-link-rate", m.config.planar_link_rate);
    pp.get("link-radius", m.config.link_radius);

    int min_wanted_d = 4;
    pp.get("d-min", min_wanted_d);

    auto samples = m.run(s);
    for (const auto& t : samples) {
        TannerGraph gr = to_tanner_graph(t);
        const uint64_t n = gr.get_vertices_by_type( tanner::vertex_t::type::data ).size();
        const uint64_t k = n - gr.get_checks().size();

        int dmin = 1000;
        std::map<int, int> dmap;
        for (int obs = 0; obs < k; obs++) {
            const int d = gr.compute_code_distance(true, obs);
            dmin = std::min(d, dmin);
            dmap[d]++;
        }
        if (dmin < min_wanted_d) continue;

        std::cout << "( rank = " << rank << " ) Sample: [[ " 
            << n << ", " << k << ", " << dmin << " ]]" << std::endl;
        std::cout << "\tDistances:\n";
        for (const auto& [d,cnt] : dmap) std::cout << "\t" << d << ": " << cnt << "\n";
    }
    MPI_Finalize();
    return 0;
}
