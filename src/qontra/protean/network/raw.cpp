/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "qontra/protean/network.h"

#include <vtils/utility.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

using namespace vtils;

RawNetwork::RawNetwork(TannerGraph& tgr) 
    :Graph(),
    v_tanner_raw_map(),
    tanner_graph(tgr),
    enable_memoization(false),
    proxy_memo_map(),
    support_memo_map(),
    same_support_memo_map(),
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

std::vector<sptr<raw_vertex_t>>&
RawNetwork::get_proxy_walk_path(sptr<raw_vertex_t> src, sptr<raw_vertex_t> dst) {

    if (!enable_memoization ||
            (!proxy_memo_map.count(src) || !proxy_memo_map[src].count(dst)))
    {
        if (contains(src, dst)) {
            tlm_put(proxy_memo_map, src, dst, {});
            tlm_put(proxy_memo_map, dst, src, {});
            goto found_proxy_path;
        }
        // We must search.
        for (sptr<raw_vertex_t> rprx : get_neighbors(src)) {
            if (rprx->qubit_type != raw_vertex_t::type::proxy) continue;
            std::vector<sptr<raw_vertex_t>> walk_path;
            if (dst == proxy_walk(src, rprx, walk_path)) {
                // We found the proxy qubit.
                tlm_put(proxy_memo_map, src, dst, walk_path);
                std::reverse(walk_path.begin(), walk_path.end());
                tlm_put(proxy_memo_map, dst, src, walk_path);
                goto found_proxy_path;
            }
        }
    }
found_proxy_path:
    return proxy_memo_map[src][dst];
}

RawNetwork::parity_support_t&
RawNetwork::get_support(sptr<raw_vertex_t> rpq) {
    if (!enable_memoization || !support_memo_map.count(rpq)) {
        parity_support_t supp;
        supp.check = rpq;
        supp.all.insert(rpq);

        // Get all data qubits in support. The flags and proxies in rpq's support are
        // any encountered when entangling the data qubits with the parity qubit.
        sptr<tanner::vertex_t> tpq = v_tanner_raw_map.at(rpq);
        for (sptr<tanner::vertex_t> tdq : tanner_graph.get_neighbors(tpq)) {
            sptr<raw_vertex_t> rdq = v_tanner_raw_map.at(tdq);
            // rdq is connected via a flag.
            if (flag_assignment_map[rpq].count(rdq)) {
                sptr<raw_vertex_t> rfq = flag_assignment_map[rpq][rdq];
                if (!contains(rfq, rdq)) {
                    insert_range(supp.proxies, get_proxy_walk_path(rfq, rdq));
                }
                if (!contains(rfq, rpq)) {
                    insert_range(supp.proxies, get_proxy_walk_path(rfq, rpq));
                }
                supp.flags.insert(rfq);
            } else if (!contains(rdq, rpq)) {
                insert_range(supp.proxies, get_proxy_walk_path(rdq, rpq));
            }
            supp.data.insert(rdq);
        }
        insert_range(supp.all, supp.data);
        insert_range(supp.all, supp.flags);
        insert_range(supp.all, supp.proxies);
        support_memo_map[rpq] = supp;
    }
    return support_memo_map[rpq];
}

sptr<net::raw_vertex_t>
RawNetwork::are_in_same_support(sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry) {
    if (!enable_memoization ||
            (!same_support_memo_map.count(rx) || !same_support_memo_map[rx].count(ry)))
    {
        // No easy way to do this besides running through all checks.
        for (sptr<raw_vertex_t> rpq : vertices) {
            if (rpq->qubit_type != raw_vertex_t::type::xparity
                    && rpq->qubit_type != raw_vertex_t::type::zparity)
            {
                continue;
            }
            parity_support_t supp = get_support(rpq);
            if (supp.all.count(rx) && supp.all.count(ry)) { 
                same_support_memo_map[rx][ry] = rpq;
                same_support_memo_map[ry][rx] = rpq;
                goto found_same_support;
            }
        }
        same_support_memo_map[rx][ry] = nullptr;
        same_support_memo_map[ry][rx] = nullptr;
    }
found_same_support:
    return same_support_memo_map[rx][ry];
}


}   // protean
}   // qontra
