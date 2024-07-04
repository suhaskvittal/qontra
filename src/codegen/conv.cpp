/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#include "codegen/conv.h"

using namespace qontra;
using namespace graph;
using namespace tanner;

namespace cct {

uptr<TannerGraph>
to_tanner_graph(uptr<TilingGraph>& gr) {
    uptr<TannerGraph> tgr = std::make_unique<TannerGraph>();

    uint64_t i = 0;
    for (sptr<shape_t> s : gr->get_vertices()) {
        auto tx = tgr->make_vertex( (i++) | VERTEX_ID_XPARITY_FLAG );
        auto tz = tgr->make_vertex( (i++) | VERTEX_ID_ZPARITY_FLAG );
        tx->qubit_type = tanner::vertex_t::type::xparity;
        tz->qubit_type = tanner::vertex_t::type::zparity;

        tgr->add_vertex(tx);
        tgr->add_vertex(tz);

        for (uint64_t q : s->qubits) {
            auto tv = tgr->get_vertex( q | VERTEX_ID_DATA_FLAG );
            if (tv == nullptr) {
                tv = tgr->make_vertex( q | VERTEX_ID_DATA_FLAG );
                tv->qubit_type = tanner::vertex_t::type::data;
                tgr->add_vertex(tv);
            }
            tgr->make_and_add_edge(tx,tv);
            tgr->make_and_add_edge(tz,tv);
        }
    }
    // Not completely done. We have to compute the logical observables.
    return tgr;
}

uptr<TannerGraph>
to_tanner_graph(const vtils::Mat2& m) {
    uptr<TannerGraph> gr = std::make_unique<TannerGraph>();
    for (uint64_t i = 0; i < m.n_cols; i++) {
        auto tv = gr->make_vertex( i | VERTEX_ID_DATA_FLAG );
        tv->qubit_type = tanner::vertex_t::type::data;
        gr->add_vertex(tv);
    }
    for (uint64_t i = 0; i < m.n_rows; i++) {
        auto tx = gr->make_vertex( (2*i) | VERTEX_ID_XPARITY_FLAG );
        auto tz = gr->make_vertex( (2*i+1) | VERTEX_ID_ZPARITY_FLAG );
        tx->qubit_type = tanner::vertex_t::type::xparity;
        tz->qubit_type = tanner::vertex_t::type::zparity;

        gr->add_vertex(tx);
        gr->add_vertex(tz);

        for (uint64_t j = 0; j < m.n_cols; j++) {
            if (m(i,j)) {
                auto tv = gr->get_vertex( j | VERTEX_ID_DATA_FLAG );
                gr->make_and_add_edge(tx,tv);
                gr->make_and_add_edge(tz,tv);
            }
        }
    }
    return gr;
}

vtils::Mat2
to_parity_matrix(uptr<TannerGraph>& gr) {
    auto data_qubits = gr->get_vertices_by_type(tanner::vertex_t::type::data);
    auto parity_checks = gr->get_vertices_by_type(tanner::vertex_t::type::xparity);
    const uint64_t nd = data_qubits.size();
    const uint64_t np = parity_checks.size();

    vtils::Mat2 m(np,nd);

    for (size_t i = 0; i < np; i++) {
        auto p = parity_checks.at(i);
        for (auto d : gr->get_neighbors(p)) {
            m.set( i, d->id & VERTEX_ID_NUMBER_MASK, true );
        }
    }
    return m;
}

}   // cct
