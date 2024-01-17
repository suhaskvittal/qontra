/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "qontra/protean/network.h"

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
    // proxy_indirection_map:
    for (auto ita = proxy_indirection_map.begin(); ita != proxy_indirection_map.end(); ) {
        // Like with flag_assignment_map, there are three possibilities for v if it exists in an entry.
        if (ita->first == v) {
            ita = proxy_indirection_map.erase(ita);
            continue;
        }
        // Goto map[*ita] (which is a map itself).
        for (auto itb = ita->second.begin(); itb != ita->second.end(); ) {
            if (itb->first == v) {  // Second key or value, so delete itb.
                itb = ita->second.erase(itb);
            } else {
                // It may be in the array.
                for (auto itc = itb->second.begin(); itc != itb->second.end(); ) {
                    if (*itc == v)  itc = itb->second.erase(itc);
                    else            itc++;
                }
                // Erase itb if itc is now empty.
                if (itb->second.empty()) itb = ita->second.erase(itb);
                else                     itb++;
            }
        }
        if (ita->second.empty()) ita = proxy_indirection_map.erase(ita);
        else                     ita++;
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
    proxy_indirection_map[prxq][q1].push_back(q2);
    proxy_indirection_map[prxq][q2].push_back(q1);
    // If q1 or q2 are proxies, we need to update their entries too.
    if (q1->qubit_type == raw_vertex_t::type::proxy) {
        update_endpoint_in_indirection_map(q1, q2, prxq);
    }
    if (q2->qubit_type == raw_vertex_t::type::proxy) {
        update_endpoint_in_indirection_map(q2, q1, prxq);
    }
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
        sptr<raw_edge_t> e = make_edge(big, rw);
        e->is_undirected = true;
        add_edge(e);
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
        proxy_proxy_merge(big, little);
    }
    delete_vertex(little);
    return little;
}

std::vector<sptr<raw_vertex_t>>
RawNetwork::proxy_walk(sptr<raw_vertex_t> from,
                        sptr<raw_vertex_t> thru,
                        std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>>& walk_res_ref) 
{
    if (!proxy_indirection_map.count(thru) || !proxy_indirection_map[thru].count(from)) {
        return {};
    }

    std::vector<sptr<raw_vertex_t>> dst_array;

    std::deque<sptr<raw_vertex_t>> bfs{thru};
    std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>> bfs_prev;
    bfs_prev[thru] = {from, thru};

    std::set<sptr<raw_vertex_t>> visited{from, thru};
    while (bfs.size()) {
        sptr<raw_vertex_t> rx = bfs.front();
        bfs.pop_front();

        if (rx->qubit_type != raw_vertex_t::type::proxy) {
            // We are done with this qubit -- it is not a proxy.
            dst_array.push_back(rx);
            walk_res_ref[rx] = std::move(bfs_prev[rx]);
            continue;
        }
        auto prev = bfs_prev[rx][bfs_prev[rx].size()-2];
        for (sptr<raw_vertex_t> ry : proxy_indirection_map[rx][prev]) {
            if (visited.count(ry)) continue;
            bfs_prev[ry] = std::vector<sptr<raw_vertex_t>>(bfs_prev[rx]);
            bfs_prev[ry].push_back(ry);
            visited.insert(ry);
            bfs.push_back(ry);
        }
    }
    return dst_array;
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

            std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>> walk_path_map;
            proxy_walk(src, rprx, walk_path_map);
            if (walk_path_map.count(dst)) {
                // We found the proxy qubit.
                auto walk_path = std::move(walk_path_map[dst]);
                tlm_put(proxy_memo_map, src, dst, walk_path);
                std::reverse(walk_path.begin(), walk_path.end());
                tlm_put(proxy_memo_map, dst, src, walk_path);
                goto found_proxy_path;
            }
        }
        std::cerr << "[ get_proxy_walk_path ] proxy path between " << print_v(src) << " and "
            << print_v(dst) << " could not be found." << std::endl;
        std::cerr << "\tadj(" << print_v(src) << "):";
        std::set<sptr<raw_vertex_t>> prx_set;
        for (sptr<raw_vertex_t> rx : get_neighbors(src)) {
            std::cerr << " " << print_v(rx);
            if (rx->qubit_type == raw_vertex_t::type::proxy) {
                std::cerr << "[";
                for (sptr<raw_vertex_t> ry : proxy_indirection_map[rx][src]) {
                    std::cerr << " " << print_v(ry);
                }
                std::cerr << " ]";
            }
        }
        std::cerr << std::endl << "\tadj(" << print_v(dst) << "):";
        for (sptr<raw_vertex_t> rx : get_neighbors(dst)) {
            std::cerr << " " << print_v(rx);
            if (rx->qubit_type == raw_vertex_t::type::proxy) {
                std::cerr << "[";
                for (sptr<raw_vertex_t> ry : proxy_indirection_map[rx][dst]) {
                    std::cerr << " " << print_v(ry);
                }
                std::cerr << " ]";
            }
        }
        std::cerr << std::endl;
        exit(1);
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
        if (rpq->qubit_type != raw_vertex_t::type::xparity
                && rpq->qubit_type != raw_vertex_t::type::zparity)
        {
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
RawNetwork::update_endpoint_in_indirection_map(
        sptr<raw_vertex_t> proxy,
        sptr<raw_vertex_t> old_endpoint,
        sptr<raw_vertex_t> new_endpoint)
{
    map_move_and_delete(proxy_indirection_map[proxy], old_endpoint, new_endpoint);
    for (auto& p : proxy_indirection_map[proxy]) {
        for (auto& x : p.second) {
            if (x == old_endpoint) x = new_endpoint;
        }
    }
}

void
RawNetwork::flag_proxy_merge(sptr<raw_vertex_t> rfq, sptr<raw_vertex_t> rprx) {
    // Get the endpoints of rprx from the proxy_indirection_map.
    for (auto& p : proxy_indirection_map[rprx]) {
        sptr<raw_vertex_t> rx = p.first;
        for (sptr<raw_vertex_t> ry : p.second) {
            // Place rx and ry in the flag_assignment_map. We also need the parity check, which
            // can be easily retrieved.
            sptr<raw_vertex_t> rpq = are_in_same_support({rfq, rx, ry});
            if (rpq == nullptr) {
                // This should not happen.
                std::cerr << "[ flag_proxy_merge ] invalid merge attempted" << std::endl;
                exit(1);
            }
            flag_assignment_map[rpq][rx] = rfq;
            flag_assignment_map[rpq][ry] = rfq;
        }
    }
    proxy_indirection_map.erase(rprx);
}

void
RawNetwork::proxy_proxy_merge(sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry) {
    // We need to join the contents of rx and ry together.
    for (auto& p : proxy_indirection_map[ry]) {
        if (rx == p.first) continue;

        if (!proxy_indirection_map[rx].count(p.first)) {
            // This is a simple move.
            proxy_indirection_map[rx][p.first] = std::move(p.second);
            continue;
        }
        // Otherwise, we need to join and handle duplicates.
        auto& arr = proxy_indirection_map[rx][p.first];
        for (sptr<raw_vertex_t> r : p.second) {
            if (std::find(arr.begin(), arr.end(), r) == arr.end()) {
                arr.push_back(r);
            }
        }
    }
    proxy_indirection_map.erase(ry);
    // We also need to remove any references to ry in rx's part of the map.
    for (auto ita = proxy_indirection_map[rx].begin(); ita != proxy_indirection_map[rx].end(); ) {
        if (ita->first == ry) {
            ita = proxy_indirection_map[rx].erase(ita);
            continue;
        }
        for (auto itb = ita->second.begin(); itb != ita->second.end(); ) {
            if (*itb == ry) itb = ita->second.erase(itb);
            else            itb++;
        }
        if (ita->second.empty()) ita = proxy_indirection_map[rx].erase(ita);
        else                     ita++;
    }
    // Now replace ry with rx everywhere else in the proxy_indirection_map
    replace_in_tracking_structures(ry, rx, 4);
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
    // flag_assignment_map:
    if (where & 2) {
        map_move_and_delete(flag_assignment_map, v, with);
        for (auto& p1 : flag_assignment_map) {
            auto& submap = p1.second;
            map_move_and_delete(submap, v, with);
            for (auto& p2 : submap) {
                if (p2.second == v) p2.second = with;
            }
        }
    }
    // proxy_indirection_map:
    if (where & 4) {
        map_move_and_delete(proxy_indirection_map, v, with);
        for (auto& p1 : proxy_indirection_map) {
            auto& submap = p1.second;
            map_move_and_delete(submap, v, with);
            for (auto& p2 : submap) {
                for (auto& x : p2.second) {
                    if (x == v) x = with;
                }
            }
        }
    }
}

}   // protean
}   // qontra
