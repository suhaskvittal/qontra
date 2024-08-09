/*
 *  author: Suhas Vittal
 *  date:   3 July 2024
 * */

#include <codegen/monte_carlo.h>
#include <codegen/convert.h>

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

    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    s = (s/size) + (rank==0)*(s%size);

    MonteCarloManager m(rank);
    m.config.r = r;
    m.config.c = c;
    auto samples = m.run(s);
    for (const auto& t : samples) {
        TannerGraph gr = to_tanner_graph(t);
        const uint64_t n = gr.get_vertices_by_type( tanner::vertex_t::type::data ).size();
        const uint64_t k = n - gr.get_checks().size();
        const int d = gr.compute_code_distance(true);

        std::cout << "( rank = " << rank << " ) Sample: [[ " 
            << n << ", " << k << ", " << d << " ]]" << std::endl;
    }
    MPI_Finalize();
    return 0;
}
