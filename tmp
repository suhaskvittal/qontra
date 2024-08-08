/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#include "codegen/driver.h"

using namespace qontra;
using namespace graph;
using namespace tanner;
using namespace vtils;

namespace cct {

void
assert_commuting(Mat2 a, Mat2 b, std::string alert) {
    for (size_t i = 0; i < a.n_rows; i++) {
        Mat2 x = a.get_row(i);
        for (size_t j = 0; j < b.n_rows; j++) {
            Mat2 y = b.get_row(j);
            auto z = x & y;
            if (z.popcount(0) & 1) {
                std::cerr << "[ " <<  alert << " ] Anticommutation detected between rows "
                    << i << " and " << j << std::endl;
                exit(1);
            }
        }
    }
}

uptr<TannerGraph>
get_sample(uint64_t max_qubits, tiling_config_t conf, int seed) {
    uptr<TilingGraph> g1 = make_random_tiling(max_qubits, conf);
    // Convert this to a parity matrix.
    uptr<TannerGraph> g2 = to_tanner_graph(g1);
    Mat2 m1 = to_parity_matrix(g2);
    std::cout << m1;
    // Check anticommutation.
    assert_commuting(m1, m1, "m1");
    // Compute the basis vectors of m. These are the stabilizer generators.
    auto v1 = get_basis_vectors(m1.trr());
    // Set the rows of m2 to the contents of v1.
    Mat2 m2 = Mat2::from_rows(v1);
    uptr<TannerGraph> g3 = to_tanner_graph(m2);
    // Now compute the logical operators.
    auto v2 = get_null_basis_vectors(m2);
    assert_commuting(m1, Mat2::from_rows(v2), "m1 vs null");

    std::cout << Mat2::from_rows(v2);

    Mat2 m3(v1.size()+v2.size(), m1.n_cols);
    m3.copy_rows_from(0, v1);
    m3.copy_rows_from(v1.size(), v2);
    // Now, we identify the logical operators by examining the pivots of m3.
    std::vector<size_t> p;
    _rref(m3.trr(), p);
    for (size_t _i : p) {
        if (_i < v1.size()) continue;
        size_t i = _i - v1.size();
        Mat2 x = v2.at(i);
        // Check if x is in the row space of m2.
        for (int b = 0; b <= 1; b++) {
            auto& obs_list = g3->get_obs_ref(b);
            std::vector<sptr<vertex_t>> obs;
            for (int j = 0; j < x.n_cols; j++) {
                if (x(0,j)) {
                    obs.push_back( g3->get_vertex( j | VERTEX_ID_DATA_FLAG ) );
                }
            }
            obs_list.push_back(obs);
        }
    }
    return g3;
}

}   // cct
