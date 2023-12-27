/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "protean/network.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

RawNetwork::Network(TannerGraph& tgr) 
    :Graph(),
    tanner_to_raw(),
    raw_to_tanner(),
    tanner_graph(tgr),
    id_ctr(0)
{
    // Essentially copy over the tanner graph.
    for (sptr<tanner::vertex_t> tv : tgr.get_vertices()) {
        sptr<raw_vertex_t> rv = std::make_shared<>();
        rv->id = id_ctr++;

        if (tv->qubit_type == tanner::vertex_t::type::xparity) {
            rv->qubit_type = raw_vertex_t::type::xparity;
        } else if (tv->qubit_type == tanner::vertex_t::type::zparity) {
            rv->qubit_type = raw_vertex_t::type::zparity;
        } else {
            rv->qubit_type = raw_vertex_t::type::data;
        }
        tanner_to_raw[tv] = rv;
        raw_to_tanner[rv] = tv;
        add_vertex(tv);
    }

    for (sptr<tanner::edge_t> te : tgr.get_edges()) {
        sptr<tanner::vertex_t> tv1 = std::reinterpret_pointer_cast<tanner::vertex_t>(te->src),
                                tv2 = std::reinterpret_pointer_cast<tanner::vertex_t>(te->dst);
        auto rv1 = tanner_to_raw[tv1],
            rv2 = tanner_to_raw[tv2];
        sptr<raw_edge_t> re = std::make_shared<>();
        re->src = rv1;
        re->dst = rv2;
        re->is_undirected = true;
        add_edge(re);
    }
}

sptr<raw_vertex_t>
RawNetwork::add_flag(sptr<raw_vertex_t> dq1, sptr<raw_vertex_t> dq2, sptr<raw_vertex_t> pq)
{
    // Get edges between the data qubits and the parity qubit.
    sptr<raw_edge_t> re1 = get_edge(dq1, pq),
                        re2 = get_edge(dq2, pq);
    delete_edge(re1);
    delete_edge(re2);
    // Create flag and corresponding edges.
    sptr<raw_vertex_t> fq = std::make_shared();
    fq->id = id_ctr++;
    fq->qubit_type = raw_vertex_t::type::flag;

    sptr<raw_edge_t> rfe1 = std::make_shared(),
                       rfe2 = std::make_shared(),
                       rfe3 = std::make_shared();
    rfe1->src = dq1;
    rfe1->dst = fq;
    rfe1->is_undirected = true;
    add_edge(rfe1);

    rfe2->src = dq2;
    rfe2->dst = fq;
    rfe2->is_undirected = true;
    add_edge(rfe2);

    rfe3->src = fq;
    rfe3->dst = pq;
    rfe3->is_undirected = true;
    add_edge(rfe3);

    // Update tracking structures and return.
    flag_ownership_map[pq].push_back(fq);
    tlm_put(flag_assignment_map, pq, dq1, fq);
    tlm_put(flag_assignment_map, pq, dq2, fq);
    return fq;
}

sptr<raw_vertex_t>
RawNetwork::add_proxy(sptr<raw_edge_t> e) {
    sptr<raw_vertex_t> q1 = std::reinterpret_pointer_cast(e->src),
                        q2 = std::reinterpret_pointer_cast(e->dst);
    delete_edge(e);
    // Create proxy and corresponding edges.
    sptr<raw_vertex_t> prxq = std::make_shared();
    prxq->id = id_ctr++;
    prxq->qubit_type = raw_vertex_t::type::proxy;

    sptr<raw_edge_t> e1 = std::make_shared(),
                        e2 = std::make_shared();
    e1->src = q1;
    e1->dst = prxq;
    e1->is_undirected = true;
    add_edge(e1);

    e2->src = prxq;
    e2->dst = q2;
    e2->is_undirected = true;
    add_edge(e2);

    // Update tracking structures and return.
    tlm_put(proxy_indirection_map, prxq, q1, q2);
    tlm_put(proxy_indirection_map, prxq, q2, q1);
    return prxq;
}

}   // protean
}   // qontra
