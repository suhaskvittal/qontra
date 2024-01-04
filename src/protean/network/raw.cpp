/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "protean/network.h"

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

RawNetwork::RawNetwork(TannerGraph& tgr) 
    :Graph(),
    v_tanner_raw_map(),
    tanner_graph(tgr),
    id_ctr(0)
{
    // Essentially copy over the tanner graph.
    for (sptr<tanner::vertex_t> tv : tgr.get_vertices()) {
        sptr<raw_vertex_t> rv = make_vertex();
        if (tv->qubit_type == tanner::vertex_t::type::xparity) {
            rv->qubit_type = raw_vertex_t::type::xparity;
        } else if (tv->qubit_type == tanner::vertex_t::type::zparity) {
            rv->qubit_type = raw_vertex_t::type::zparity;
        } else {
            rv->qubit_type = raw_vertex_t::type::data;
        }
        v_tanner_raw_map.put(tv, rv);
        add_vertex(rv);
    }

    for (sptr<tanner::edge_t> te : tgr.get_edges()) {
        sptr<tanner::vertex_t> tv1 = tgr.get_source(te),
                                tv2 = tgr.get_target(te);
        auto rv1 = v_tanner_raw_map.at(tv1),
            rv2 = v_tanner_raw_map.at(tv2);
        sptr<raw_edge_t> re = make_edge(rv1, rv2);
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
    sptr<raw_vertex_t> fq = make_vertex();
    fq->qubit_type = raw_vertex_t::type::flag;
    add_vertex(fq);

    sptr<raw_edge_t> rfe1 = make_edge(dq1, fq),
                       rfe2 = make_edge(dq2, fq),
                       rfe3 = make_edge(fq, pq);
    add_edge(rfe1);
    add_edge(rfe2);
    add_edge(rfe3);

    // Update tracking structures and return.
    flag_ownership_map[pq].push_back(fq);
    tlm_put(flag_assignment_map, pq, dq1, fq);
    tlm_put(flag_assignment_map, pq, dq2, fq);
    return fq;
}

sptr<raw_vertex_t>
RawNetwork::add_proxy(sptr<raw_edge_t> e) {
    sptr<raw_vertex_t> q1 = get_source(e),
                        q2 = get_target(e);
    delete_edge(e);
    // Create proxy and corresponding edges.
    sptr<raw_vertex_t> prxq = make_vertex();
    prxq->qubit_type = raw_vertex_t::type::proxy;
    add_vertex(prxq);

    sptr<raw_edge_t> e1 = make_edge(q1, prxq),
                        e2 = make_edge(prxq, q2);
    add_edge(e1);
    add_edge(e2);
    // Update tracking structures and return.
    tlm_put(proxy_indirection_map, prxq, q1, q2);
    tlm_put(proxy_indirection_map, prxq, q2, q1);
    return prxq;
}

}   // protean
}   // qontra
