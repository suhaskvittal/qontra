/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "qontra/protean/network.h"
#include "qontra/graph/algorithms/distance.h"

#include <vtils/utility.h>

#include <algorithm>
#include <deque>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

using namespace vtils;

inline int
get_merge_priority(sptr<raw_vertex_t> v) {
    if (v->qubit_type == raw_vertex_t::type::xparity
            || v->qubit_type == raw_vertex_t::type::zparity) 
    {
        return 2;
    } else if (v->qubit_type == raw_vertex_t::type::flag) {
        return 1;
    } else {
        return 0;
    }
}

// Moves the entry of v to that of w and erases v.
template <class T, class MAP> inline void
map_move_and_delete(MAP& m, T v, T w) {
    if (!m.count(v)) return;
    m[w] = std::move(m[v]);
    m.erase(v);
}

RawNetwork::RawNetwork(TannerGraph* tgr) 
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
    for (sptr<tanner::vertex_t> tv : tgr->get_vertices()) {
        sptr<raw_vertex_t> rv = make_and_add_vertex();
        if (tv->qubit_type == tanner::vertex_t::type::xparity) {
            rv->qubit_type = raw_vertex_t::type::xparity;
        } else if (tv->qubit_type == tanner::vertex_t::type::zparity) {
            rv->qubit_type = raw_vertex_t::type::zparity;
        } else {
            rv->qubit_type = raw_vertex_t::type::data;
        }
        v_tanner_raw_map.put(tv, rv);
    }

    for (sptr<tanner::edge_t> te : tgr->get_edges()) {
        sptr<tanner::vertex_t> tv1 = te->get_source<tanner::vertex_t>(),
                                tv2 = te->get_target<tanner::vertex_t>();
        auto rv1 = v_tanner_raw_map.at(tv1),
            rv2 = v_tanner_raw_map.at(tv2);
        sptr<raw_edge_t> re = make_and_add_edge(rv1, rv2);
    }
}

void
RawNetwork::delete_vertex(sptr<raw_vertex_t> v) {
    Graph::delete_vertex(v);
    // Remove v from any data structures.
    if (v_tanner_raw_map.count(v)) {
        v_tanner_raw_map.erase(v);
    }
    // flag_ownership_map:
    for (auto ita = flag_ownership_map.begin(); ita != flag_ownership_map.end(); ) {
        // If it is a key, remove the entry.
        if (ita->first == v) {
            ita = flag_ownership_map.erase(ita);
            continue;
        }
        // Otherwise, remove it from the array of flags (if it exists).
        for (auto itb = ita->second.begin(); itb != ita->second.end(); ) {
            if (*itb == v)  itb = ita->second.erase(itb);
            else            itb++;
        }
        // Delete ita if the array of flags becomes empty.
        if (ita->second.empty()) ita = flag_ownership_map.erase(ita);
        else                     ita++;
    }
    // flag_assignment_map:
    for (auto ita = flag_assignment_map.begin(); ita != flag_assignment_map.end(); ) {
        // Three possibilities if v is in an entry: (1) first key, (2) second key, or (3) value.
        if (ita->first == v) {  // It is the first key, so delete ita.
            ita = flag_assignment_map.erase(ita);
            continue;
        }
        // Goto map[*ita] (which is a map itself).
        for (auto itb = ita->second.begin(); itb != ita->second.end(); ) {
            if (itb->first == v || itb->second == v) {  // Second key or value, so delete itb.
                itb = ita->second.erase(itb);
            } else {
                itb++;
            }
        }
        if (ita->second.empty()) ita = flag_assignment_map.erase(ita);
        else                     ita++;
    }
    // flag_support_map.
    for (auto ita = flag_support_map.begin(); ita != flag_support_map.end(); ) {
        if (ita->first == v) {
            ita = flag_support_map.erase(ita);
            continue;
        }
        // map[*ita], which is itself a map.
        for (auto itb = ita->second.begin(); itb != ita->second.end(); ) {
            if (itb->first == v) {
                itb = ita->second.erase(itb);
            } else {
                // itb->second is a set, so we can try to erase v from the set.
                itb->second.erase(v);
                if (itb->second.empty()) itb = ita->second.erase(itb);
                else itb++;
            }
        }
        if (ita->second.empty()) ita = flag_support_map.erase(ita);
        else ita++;
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
    sptr<raw_vertex_t> fq = make_and_add_vertex();
    fq->qubit_type = raw_vertex_t::type::flag;

    make_and_add_edge(dq1, fq);
    make_and_add_edge(dq2, fq);
    make_and_add_edge(fq, pq);

    // Update tracking structures and return.
    flag_ownership_map[pq].push_back(fq);
    tlm_put(flag_assignment_map, pq, dq1, fq);
    tlm_put(flag_assignment_map, pq, dq2, fq);
    tlm_put(flag_support_map, pq, fq, std::set<sptr<raw_vertex_t>>{dq1, dq2});
    if (pq->qubit_type == raw_vertex_t::type::xparity) {
        x_flag_set.insert(fq);
    }
    return fq;
}

sptr<raw_vertex_t>
RawNetwork::add_proxy(sptr<raw_edge_t> e) {
    sptr<raw_vertex_t> q1 = e->get_source<raw_vertex_t>(),
                        q2 = e->get_target<raw_vertex_t>();
    delete_edge(e);
    // Create proxy and corresponding edges.
    sptr<raw_vertex_t> prxq = make_and_add_vertex();
    prxq->qubit_type = raw_vertex_t::type::proxy;

    make_and_add_edge(q1, prxq);
    make_and_add_edge(prxq, q2);
    return prxq;
}

sptr<raw_vertex_t>
RawNetwork::merge(sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry) {
    int rx_prio = get_merge_priority(rx),
        ry_prio = get_merge_priority(ry);
    // The lower priority (less important) role will be merged into the higher priority role.
    sptr<raw_vertex_t> big, little;
    if (rx_prio >= ry_prio) { big = rx; little = ry; }
    else                    { big = ry; little = rx; }
    // Perform the merge.
    for (sptr<raw_vertex_t> rw : get_neighbors(little)) {
        if (big == rw) continue;
        if (contains(big, rw)) continue;
        make_and_add_edge(big, rw);
    }
    // Depending on what big and little are, we need to update the
    // tracking data structures.
    if (big->qubit_type == raw_vertex_t::type::xparity || big->qubit_type == raw_vertex_t::type::zparity) {
        replace_in_tracking_structures(little, big, 4);
    } else if (big->qubit_type == raw_vertex_t::type::flag) {
        if (little->qubit_type == raw_vertex_t::type::flag) {
            replace_in_tracking_structures(little, big);
        } else {
            flag_proxy_merge(big, little);
        }
    } else {
        replace_in_tracking_structures(little, big);
    }
    delete_vertex(little);
    return little;
}


std::vector<sptr<raw_vertex_t>>&
RawNetwork::get_proxy_walk_path(sptr<raw_vertex_t> src, sptr<raw_vertex_t> dst) {

    if (!enable_memoization ||
            (!proxy_memo_map.count(src) || !proxy_memo_map[src].count(dst)))
    {
        if (contains(src, dst)) {
            tlm_put(proxy_memo_map, src, dst, {src, dst});
            tlm_put(proxy_memo_map, dst, src, {dst, src});
            goto found_proxy_path;
        }
        // Compute Dijkstra's from src.
        std::map<sptr<raw_vertex_t>, fp_t> dist;
        std::map<sptr<raw_vertex_t>, sptr<raw_vertex_t>> pred;
        dijkstra(this, src, dist, pred,
            [&] (sptr<raw_edge_t> e)
            {
                sptr<raw_vertex_t> x = e->get_source<raw_vertex_t>(),
                                    y = e->get_target<raw_vertex_t>();
                if (x->qubit_type == raw_vertex_t::type::data
                    || y->qubit_type == raw_vertex_t::type::data)
                {
                    return 100.0;
                }
                return 1.0;
            });
        // If enable_memoization = true, then compute distance to all vertices. Otherwise,
        // only compute to dst.
        std::vector<sptr<raw_vertex_t>> compute_distance_to{dst};
        if (enable_memoization) compute_distance_to = vertices;
        // Now, use pred to compute paths.
        for (sptr<raw_vertex_t> rx : compute_distance_to) {
            sptr<raw_vertex_t> curr = rx;
            std::vector<sptr<raw_vertex_t>> path;
            while (curr != src) {
                path.push_back(curr);
                curr = pred.at(curr);
            }
            path.push_back(src);
            // Well path is currently oriented from rx to src.
            tlm_put(proxy_memo_map, rx, src, path);
            std::reverse(path.begin(), path.end());
            tlm_put(proxy_memo_map, src, rx, path);
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
        for (sptr<tanner::vertex_t> tdq : tanner_graph->get_neighbors(tpq)) {
            sptr<raw_vertex_t> rdq = v_tanner_raw_map.at(tdq);
            // rdq is connected via a flag.
            if (flag_assignment_map[rpq].count(rdq)) {
                sptr<raw_vertex_t> rfq = flag_assignment_map[rpq][rdq];
                if (!contains(rfq, rdq)) {
                    insert_range(supp.proxies, get_proxy_walk_path(rfq, rdq), 1, -1);
                }
                if (!contains(rfq, rpq)) {
                    insert_range(supp.proxies, get_proxy_walk_path(rfq, rpq), 1, -1);
                }
                supp.flags.insert(rfq);
            } else if (!contains(rdq, rpq)) {
                insert_range(supp.proxies, get_proxy_walk_path(rdq, rpq), 1, -1);
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

sptr<raw_vertex_t>
RawNetwork::are_in_same_support(std::initializer_list<sptr<raw_vertex_t>> vertex_list) {
    for (sptr<raw_vertex_t> rpq : vertices) {
        if (!rpq->is_check()) {
            continue;
        }
        parity_support_t supp = get_support(rpq);
        bool all_in_supp = true;
        for (sptr<raw_vertex_t> rv : vertex_list) {
            all_in_supp &= supp.all.count(rv);
        }
        if (all_in_supp) return rpq;
    }
    return nullptr;
}

void
RawNetwork::flag_proxy_merge(sptr<raw_vertex_t> rfq, sptr<raw_vertex_t> rprx) {
    // Get the neighbors of rprx that are not rfq
    for (sptr<raw_vertex_t> rx : get_neighbors(rprx)) {
        if (rx->qubit_type == raw_vertex_t::type::data) {
            // Then just add to the flag_assignment_map
            sptr<raw_vertex_t> rpq = are_in_same_support({rfq, rx});
            if (rpq == nullptr) {
                // This should not happen.
                std::cerr << "[ flag_proxy_merge ] invalid merge attempted" << std::endl;
                exit(1);
            }
            flag_assignment_map[rpq][rx] = rfq;
            flag_support_map[rpq][rfq].insert(rx);
        }
    }
}

void
RawNetwork::replace_in_tracking_structures(sptr<raw_vertex_t> v, sptr<raw_vertex_t> with, uint8_t where) {
    // flag_ownership_map:
    if (where & 1) {
        map_move_and_delete(flag_ownership_map, v, with);
        for (auto& p1 : flag_ownership_map) {
            for (auto& x : p1.second) {
                if (x == v) x = with;
            }
        }
    }
    // flag_assignment_map and flag_support_map:
    if (where & 2) {
        map_move_and_delete(flag_assignment_map, v, with);
        for (auto& [_tmp, submap] : flag_assignment_map) {
            map_move_and_delete(submap, v, with);
            for (auto& [_tmp, x] : submap) {
                if (x == v) x = with;
            }
        }
        map_move_and_delete(flag_support_map, v, with);
        for (auto& [_tmp, submap]: flag_support_map) {
            map_move_and_delete(submap, v, with);
            for (auto& [_tmp, x] : submap) {
                if (x.count(v)) {
                    x.erase(v);
                    x.insert(with);
                }
            }
        }
    }
}

}   // protean
}   // qontra
