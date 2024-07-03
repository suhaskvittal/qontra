/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#include "codegen/driver.h"

using namespace qontra;
using namespace graph;

namespace cct {

uptr<TannerGraph>
get_sample(uint64_t max_qubits, tiling_config_t conf, int seed) {
    uptr<TilingGraph> g1 = make_random_tiling(max_qubits, conf);
    // Convert this to a parity matrix.
    uptr<TannerGraph> g2 = to_tanner_graph(g1);
    vtils::Mat2 m1 = to_parity_matrix(g2);
    // Compute the basis vectors of m. These are the stabilizer generators.
    auto v1 = vtils::get_basis_vectors(m.trr());
    // Set the rows of m2 to the contents of v1.
    vtils::Mat2 m2(v1.size(), m1.n_cols);
    for (size_t i = 0; i < v1.size(); i++) {
        for (size_t j = 0; j < m1.n_cols; j++) {
            m2.set(i, j, v1[i](0,j));
        }
    }
    uptr<TannerGraph> g3 = to_tanner_graph(m2);
    // Now compute the logical operators.
    auto v2 = vtils::get_null_basis_vectors(m);
    vtils::Mat2 m3(v1.size()+v2.size(), m1.n_cols);
    for (size_t i = 0; i < v1.size(); i++) {
        for (size_t j = 0; j < m1.n_cols; j++) {
            m3.set(i, j, v1[i](0,j));
        }
    }
    for (size_t i = 0; i < v2.size(); i++) {
        for (size_t j = 0; j < m1.n_cols; j++) {
            m3.set(i+v1.size(), j, v1[i](0,j));
        }
    }
    // Now, we identify the logical operators by examining the pivots of m3.
    std::vector<size_t> p;
    vtils::_rref(m3, p);
    for (size_t _i : p) {
        if (_i < v1.size()) continue;
        size_t i = _i - v1.size();
        Mat2 x = v2.at(i);
        for (int b = 0; b <= 1; b++) {
            auto& obs_list = g3->get_obs_ref(b);
            std::vector<sptr<tanner::vertex_t>> obs;
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
