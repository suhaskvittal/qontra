/*
 *  author: Suhas Vittal
 *  date:   9 August 2024
 * */
#define DEBUG

#include "codegen/convert.h"

#include <vtils/set_algebra.h>

using namespace qontra;
using namespace graph;
using namespace tanner;

using namespace vtils;

namespace cgen {

ctri_t
make_ctri(sptr<check_t> x, int sx, sptr<check_t> y, int sy, sptr<check_t> z, int sz) {
    std::array<sptr<check_t>, 3> xyz{x,y,z};
    std::map<sptr<check_t>, int> sides{
        std::make_pair(x,sx),
        std::make_pair(y,sy),
        std::make_pair(z,sz)
    };
    std::sort(xyz.begin(), xyz.end());
    return std::make_tuple( xyz[0], sides[xyz[0]],
                            xyz[1], sides[xyz[1]],
                            xyz[2], sides[xyz[2]] );
}

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
                std::cerr << "\t" << x << "\n\t" << y << "\n\t" << z << "\n";
                exit(1);
            }
        }
    }
}

TannerGraph
to_tanner_graph(const Tiling& t) {
    Mat2 m1 = to_matrix(t);

    assert_commuting(m1, m1, "m1");

    auto b1 = get_basis_vectors(m1.trr());  // stabilizer generators.
    Mat2 m2 = Mat2::from_rows(b1);
    auto b2 = get_null_basis_vectors(m2);  // logical operators.

    assert_commuting(m1, Mat2::from_rows(b2), "m1 vs null");

    Mat2 m3(b1.size()+b2.size(), m1.n_cols);
    m3.copy_rows_from(0, b1);
    m3.copy_rows_from(b1.size(), b2);
    std::vector<size_t> p;
    _rref(m3.trr(), p);

    TannerGraph gr = to_tanner_graph(m2);
    for (size_t _i : p) {
        if (_i < b1.size()) continue;
        size_t i = _i - b1.size();
        Mat2 x = b2.at(i);
        // Check if x is in the row space of m2.
        for (int b = 0; b <= 1; b++) {
            auto& obs_list = gr.get_obs_ref(b);
            std::vector<sptr<vertex_t>> obs;
            for (int j = 0; j < x.n_cols; j++) {
                if (x(0,j)) {
                    obs.push_back( gr.get_vertex( j | VERTEX_ID_DATA_FLAG ) );
                }
            }
            obs_list.push_back(obs);
        }
    }
    return gr;
}

Mat2
to_matrix(const Tiling& t) {
    const auto& checks = t.get_all_checks_ref();

    std::set<ctri_t> qubits;
    for (sptr<check_t> c_cen : t.get_checks_ref(0)) {
//      std::cout << print_check(c_cen) << " qubits:\n";
        for (int s = 0; s < c_cen->size(); s++) {
            sptr<check_t> c_left = c_cen->get(s),
                          c_right = c_cen->get(s+1);
            int s_left = (c_left==nullptr) ? -1 : c_left->get_const_side_of(c_cen, s),
                s_right = (c_right==nullptr) ? -1 : c_right->get_const_side_of(c_cen, s+1);

            ctri_t t = make_ctri(c_left, s_left, c_cen, s, c_right, s_right);
//          std::cout << "\t" << print_ctri(t) << " " << qubits.count(t) << "\n";
            qubits.insert(t);

        }
    }

    const int n = qubits.size();
    const int m = checks.size();
//  std::cout << "n, m = " << n << ", " << m << "\n";
    Mat2 mat(m, n);
    // Make enumeration map for checks.
    std::map<sptr<check_t>, int> check_enum_map;
    for (int i = 0; i < checks.size(); i++) check_enum_map[checks.at(i)] = i;
    int k = 0;
    for (const auto& q : qubits) {
        const auto& [x,sx,y,sy,z,sz] = q;
//      std::cout << "Qubit " << k << " = " << print_ctri(q);
        for (sptr<check_t> c : {x,y,z}) {
            if (c != nullptr) {
                mat.set(check_enum_map[c], k, true);
            }
        }
//      std::cout << "\n";
        k++;
    }
    return mat;
}

TannerGraph
to_tanner_graph(const Mat2& mat) {
    TannerGraph gr;
    for (uint64_t i = 0; i < mat.n_cols; i++) {
        auto tv = gr.make_vertex( i | VERTEX_ID_DATA_FLAG );
        tv->qubit_type = tanner::vertex_t::type::data;
        gr.add_vertex(tv);
    }
    for (uint64_t i = 0; i < mat.n_rows; i++) {
        auto tx = gr.make_vertex( (2*i) | VERTEX_ID_XPARITY_FLAG );
        auto tz = gr.make_vertex( (2*i+1) | VERTEX_ID_ZPARITY_FLAG );
        tx->qubit_type = tanner::vertex_t::type::xparity;
        tz->qubit_type = tanner::vertex_t::type::zparity;

        gr.add_vertex(tx);
        gr.add_vertex(tz);

        for (uint64_t j = 0; j < mat.n_cols; j++) {
            if (mat(i,j)) {
                auto tv = gr.get_vertex( j | VERTEX_ID_DATA_FLAG );
                gr.make_and_add_edge(tx,tv);
                gr.make_and_add_edge(tz,tv);
            }
        }
    }
    return gr;
}

}   // cgen
