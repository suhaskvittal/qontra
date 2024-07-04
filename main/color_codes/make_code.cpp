/*
 *  author: Suhas Vittal
 *  date:   3 July 2024
 * */

#include <codegen/driver.h>

#include <vtils/cmd_parse.h>

#include <string>

#include <stdlib.h>

using namespace qontra;
using namespace graph;
using namespace cct;
using namespace vtils;

int main(int argc, char* argv[]) {
    uint64_t nmax = atoll(argv[1]);

    tiling_config_t conf;
    uptr<TannerGraph> gr = get_sample(nmax, conf, 0);
    // Compute code params.
    const uint64_t n = gr->get_vertices_by_type( tanner::vertex_t::type::data ).size();
    const uint64_t k = n - gr->get_checks().size();
    const int d = gr->compute_code_distance(true);

    std::cout << "[[n, k, d]] = [[" << n << ", " << k << ", " << d << "]]\n";
}
