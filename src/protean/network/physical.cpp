/* author: Suhas Vittal
 *  date:   27 December 2023
 * */


#include "protean/network.h"

#include <defs/set_algebra.h>

#include <lemon/matching.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

template <class T> inline std::pair<T, T>
make_ordered_pair(T x, T y) {
    if (x < y)  return std::make_pair(x, y);
    else        return std::make_pair(y, x);
}

template <class CONTAINER, class T> inline std::set<T>
make_set(CONTAINER<T> x) {
    return std::set<T>(x.begin(), x.end());
}

PhysicalNetwork(TannerGraph& tgr)
    :Graph(),
    raw_connection_network(tgr),
    role_to_phys(),
    // Planarity tracking:
    planar_repr(),
    v_phys_lemon_map(),
    e_phys_lemon_map(),
    // Other tracking:
    occupied_tsvs(),
    bulk_degree_map(),
    // Other variables:
    id_ctr(0),
{
    // Create a corresponding physical qubit for every vertex in raw_connection_network.
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        sptr<phys_vertex_t> pv = std::make_shared<>();
        pv->id = id_ctr++;
        pv->push_back_role(rv);
        role_to_phys[rv] = pv;
        add_vertex(pv);
    }
    // And the corresponding edges as well.
    for (sptr<raw_edge_t> re : raw_connection_network.get_edges()) {
        sptr<raw_vertex_t> rsrc = std::reinterpret_pointer_cast<raw_vertex_t>(re->src),
                            rdst = std::reinterpret_pointer_cast<raw_vertex_t>(re->dst);
        sptr<phys_edge_t> pe = std::make_shared<>();
        pe->src = role_to_phys[rsrc];
        pe->dst = role_to_phys[rdst];
        pe->is_undirected = true;
        add_edge(pe);
    }
}

void
PhysicalNetwork::join_qubits_with_identical_support() {
    // These are vertices that are consumed by another.
    std::set<sptr<phys_vertex_t>> deleted_vertices;

    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<phys_vertex_t> pv = vertices[i];
        if (deleted_vertices.count(pv)) continue;

        std::set<sptr<phys_vertex_t>> pv_adj = make_set(get_neighbors(pv));
        for (size_t j = i+1; j < vertices.size(); j++) {
            sptr<phys_vertex_t> pw = vertices[j];
            if (deleted_vertices.count(pw)) continue;

            size_t common_neighbors = 0;
            for (sptr<phys_vertex_t> x : pv_adj) {
                common_neighbors += contains(pw, x);
            }
            if (common_neighbors == pv_adj.size()) {
                pv->consume(pw);    // Now, pv contains all the roles of pw.
                deleted_vertices.insert(pw);
            }
        }
    }
    // Delete the consumed vertices
    for (sptr<phys_vertex_t> x : deleted_vertices) delete_vertex(x);
    reallocate_edges();
}

void
PhysicalNetwork::make_flags() {
    using namespace lemon;
    // Algorithm:
    //  High-level: saturate each check with the maximum number of flags, and then remove
    //              any unnecessary flags via error propagation analysis.
    //  (1) For each check:
    //      (i) Assign flags in a way that maximizes pairs of qubits with many common
    //          checks.
    //      (ii) If any pair of qubits were assigned for a prior check, pre-assign them.
    //  (2) After assigning flags, identify any flags that do not protect against weight-2
    //      errors. Remove these flags.
    //  (3) For the remaining flags (which are raw vertices), create corresponding physical
    //      vertices and update the connectivity.
    TannerGraph& tanner_graph = raw_connection_network.tanner_graph;

    template <class V> using v_pair_t=std::pair<sptr<V>, sptr<V>>;

    // The entries should be ordered pairs to avoid double counting.
    typedef std::set<v_pair_t<raw_vertex_t>>    flag_pair_set_t;

    std::map<raw_vertex_t, flag_pair_set_t> proposed_flag_pair_map;
    flag_pair_set_t all_proposed_flag_pairs;

    for (sptr<tanner::vertex_t> check : tanner_graph.get_checks()) {
        std::vector<sptr<tanner::vertex_t>> _support = tanner_graph.get_neighbors(check);
        // Get raw vertices.
        sptr<raw_vertex_t> rpq = raw_connection_network.v_tanner_raw_map.at(check);
        std::vector<sptr<raw_vertex_t>> support;
        for (sptr<tanner::vertex_t> tx : _support) {
            sptr<raw_vertex_t> rx = raw_connection_network.v_tanner_raw_map.at(tx);
            support.push_back(rx);
        }
        // Now, compute the flag setup.
        ListGraph gr;   // This is for max-weight matching.
        ListGraph::EdgeMap<int> w_map;

        BijectiveMap<sptr<raw_vertex_t>, ListGraph::Node> data_lemon_map;
        // First, create a node for each vertex in the support.
        for (sptr<raw_vertex_t> rv : support) data_lemon_map.put(rv, gr.addNode());

        flag_pair_set_t flag_pairs;

        std::set<raw_vertex_t> visited;
        for (size_t i = 0; i < support.size(); i++) {
            sptr<raw_vertex_t> rv = support[i];
            if (visited.count(rv)) continue;

            auto vnode = data_lemon_map.at(rv);
            // We will the physical qubit (specifically its neighbors) each tick of the
            // inner loop.
            //
            // Note that we only care about neighbors that are parity qubits.
            sptr<phys_vertex_t> pv = role_to_phys[rv];
            std::set<sptr<phys_vertex_t>> p_neighbors;
            for (sptr<phys_vertex_t> x : get_neighbors(pv)) {
                if (x->has_role_of_type(raw_vertex_t::type::xparity)
                    || x->has_role_of_type(raw_vertex_t::type::zparity)
                {
                    p_neighbors.insert(x);
                }
            }
            for (size_t j = i+1; j < support.size(); j++) {
                sptr<raw_vertex_t> rw = support[j];
                if (visited.count(rw)) continue;
            
                if (proposed_flag_pair_map.count(make_ordered_pair(rv, rw))) {
                    // Then, rv and rw are already assigned. Move on.
                    visited.insert(rv);
                    visited.insert(rw);
                    // Erase the lemon vertices so no vertex can match to rv or rw.
                    gr.erase(data_lemon_map.at(rv));
                    gr.erase(data_lemon_map.at(rw));
                    // Add this pair to the flag pairs.
                    flag_pairs.insert(make_ordered_pair(rv, rw));
                    continue;
                }

                auto wnode = data_lemon_map.at(rw);
                auto edge = gr.addEdge(vnode, wnode);
                // Get corresponding physical vertices, get number of common neighbors,
                // and set the weight of edge to that.
                sptr<phys_vertex_t> pw = role_to_phys[rw];
                size_t common_neighbors = 0;
                for (sptr<phys_vertex_t> x : p_neighbors) {
                    common_neighbors += contains(pw, x);
                }
                w_map[edge] = static_cast<int>(common_neighbors);
            }
        }
        // Now, perform Max Weight Matching.
        MaxWeightPerfectMatching maxwpm(gr, w_map);
        maxwpm.run();
        for (sptr<raw_vertex_t> rv : support) {
            if (visited.count(rv)) continue;

            auto vnode = data_lemon_map.at(rv);
            auto wnode = maxwpm.mate(vnode);
            sptr<raw_vertex_t> rw = data_lemon_map.at(wnode);
            
            flag_pairs.insert(make_ordered_pair(rv, rw));
            visited.insert(rv);
            visited.insert(rw);
        }
        // Update the flag proposal structures.
        proposed_flag_pair_map[rpq] = flag_pairs;
        all_proposed_flag_pairs += flag_pairs;
    }
    // Now, we must run an error-propagation analysis to identify any useless flags.
    // 
    // Separate flags from those that protect against X and Z errors.
    flag_pair_set_t x_flags, z_flags;
    for (auto& pair : proposed_flag_pair_map) {
        sptr<raw_vertex_t> rpq = pair.first;
        flag_pair_set_t flags = pair.second;
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            x_flags += flags;
        } else {
            z_flags += flags;
        }
    }
    std::array<flag_pair_set_t, 2> flag_pair_sets{x_flags, z_flags};
    for (size_t i = 0; i < 2; i++) {
        flag_pair_set_t& flags = flag_pair_sets[i];
        stim::simd_bits<SIMD_WIDTH> indicator_bits = do_flags_protect_weight_two_error(flags, i==0);
        // Remove any flags whose indicator bit is 0.
        size_t k = 0;
        for (auto it = flags.begin(); it != flags.end(); ) {
            if (indicator_bits[k++]) {
                it++;
            } else {
                all_proposed_flag_pairs -= *it;
                it = flags.erase(it);
            }
        }
    }

    // First, we will create flags for each parity qubit. This ignores the case where (q1, q2)
    // is a flag in two different parity qubits. However, note we are only interested in establishing
    // roles here. These common roles will eventually be compacted into one flag.
    //
    // We track common roles in the following structure. The values are a parity-qubit and flag-qubit pair.
    std::map<v_pair_t<raw_vertex_t>, std::set<v_pair_t<raw_vertex_t>>>
        flag_pairs_to_flag_roles;

    for (auto& pair : proposed_flag_pair_map) {
        sptr<raw_vertex_t> rpq = pair.first;
        flag_pair_set_t flags = pair.second;
        flags *= all_proposed_flag_pairs;   // Any flags that are not in the all_proposed_flag_pairs
                                            // set will be removed.
        for (auto& fp : flags) {
            sptr<raw_vertex_t> rv1 = fp.first,
                                rv2 = fp.second;
            sptr<raw_vertex_t> rfq = raw_connection_network.add_flag(rv1, rv2, rpq);
            // Update flag_pair_to_flag_roles
            flag_pairs_to_flag_roles[fp].insert(std::make_pair(rpq, rfq));
        }
    }
    // Now, we can make the physical qubits and update the connectivity.
    for (auto& p1 : flag_pairs_to_flag_roles) {
        sptr<raw_vertex_t> rv1 = p1.first.first,
                            rv2 = p1.first.second;
        std::set<v_pair_t<raw_vertex_t>> role_set = p1.second;

        sptr<phys_vertex_t> pv1 = role_to_phys[rv1],
                            pv2 = role_to_phys[rv2];

        sptr<phys_vertex_t> pfq = std::make_shared<>();
        pfq->id = id_ctr++;

        std::set<sptr<phys_vertex_t>> qubits_connected_to_pfq{pv1, pv2};
        for (auto& p2 : role_set) {
            sptr<raw_vertex_t> rpq = p2.first,
                                rfq = p2.second;
            sptr<phys_vertex_t> ppq = role_to_phys[rpq];

            // Compute cycle for role. This is the maximum of the cycles of rv1, rv2, and rfq.
            size_t cycle = std::max(pv1->cycle_role_map.at(rv1),
                                    pv2->cycle_role_map.at(rv2),
                                    ppq->cycle_role_map.at(rpq));
            pfq->add_role(rfq, cycle);
            parity_qubits_connected_to_pfq.insert(ppq);
            // Delete edge between pv1/pv2 and ppq (if it exists).
            if (contains(pv1, ppq)) {
                delete_edge(get_edge(pv1, ppq));
            } 
            if (contains(pv2, ppq)) {
                delete_edge(get_edge(pv1, ppq));
            }
        }
        // Add edges for pfq.
        for (sptr<phys_vertex_t> x : qubits_connected_to_pfq) {
            sptr<phys_edge_t> pe = std::make_shared<>();
            pe->src = pfq;
            pe->dst = x;
            pe->is_undirected = true;
            add_edge(pe);
        }
    }
}

void
PhysicalNetwork::add_connectivity_reducing_proxies() {
    bool were_any_proxies_added = false;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        const size_t dg = get_degree(pv);
        if (dg <= config.max_connectivity) continue;

        were_any_proxies_added = true;

        size_t slack_violation = dg - config.max_connectivity;
        std::vector<sptr<phys_vertex_t>> adj = get_neighbors(pv);

        sptr<phys_vertex_t> pprx = std::make_shared<>();
        pprx->id = id_ctr++;
        // Let k = slack_violation. Then the top k neighbors will
        // be serviced by the proxy.
        std::set<sptr<phys_vertex_t>> share_an_edge_with_pprx{pv};
        for (size_t i = 0; i < slack_violation; i++) {
            sptr<phys_vertex_t> px = adj[i];
            // Now, for each role in pv and px, add a proxy in the raw_connection_network.
            // Each added proxy will be a role for pprx.
            for (sptr<raw_vertex_t> rv : pv->role_set) {
                for (sptr<raw_vertex_t> rx : px->role_set) {
                    sptr<raw_edge_t> re = raw_connection_network.get_edge(rv, rx);
                    if (re == nullptr) continue;
                    sptr<raw_vertex_t> rprx = raw_connection_network.add_proxy(re);
                    // Get cycle of role (which is just the max of the cycles of rv and rx).
                    size_t cycle = std::max(pv->cycle_role_map.at(rv), px->cycle_role_map.at(rx));
                    pprx.add_role(rprx, cycle);
                }
            }
            share_an_edge_with_pprx.insert(px);
            // Delete px's edge with pv.
            sptr<phys_edge_t> pe = get_edge(pv, px);
            delete_edge(pv, px);
        }
        // Update the physical connectivity on our side.
        for (sptr<phys_vertex_t> x : share_an_edge_with_pprx) {
            sptr<phys_edge_t> e = std::make_shared<>();
            e->src = pprx;
            e->dst = x;
            e->is_undirected = true;
            add_edge(e);
        }
    }
    // If any proxies were added, we need to check the connectivity of these proxies.
    if (were_any_proxies_added) add_connectivity_reducing_proxies();
}

stim::simd_bits<SIMD_WIDTH>
PhysicalNetwork::do_flags_protect_weight_two_error(
        std::set<std::pair<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>> flag_pairs,
        bool is_x_error) 
{
    stim::simd_bits<SIMD_WIDTH> result(flag_pairs.size());
    result.invert_bits();   // All should be set to 1 to start with.
    // Algorithm: we need to check against any logical operator. This can be rather exhaustive.
    // We will use the concept of an operator tree.
    // 
    // The lowest level of the tree corresponds to the observables of the code.
    TannerGraph& tanner_graph = raw_connection_network.tanner_graph;
    auto obs_list = tanner_graph.get_obs(is_x_error);

    // We will only examine a single level of the tree at any time to avoid high memory overheads of
    // loading the entire tree.
    const size_t OPERATOR_TREE_SIZE_LIMIT = 1024L*1024L;

    // The first entry are the qubits in the operator.
    // THe second entry is a simple hash that allows us to check if two operators are equivalent.
    typedef std::tuple<<std::set<sptr<raw_vertex_t>>, stim::simd_bits<SIMD_WIDTH>> quop_t;
    std::vector<quop_t> curr_operator_tree_level;

    const size_t HASH_SIZE = obs_list.size() + tanner_graph.get_checks().size()/2;
    const size_t HASH_STABILIZER_OFFSET = obs_list.size();
    // Populate the lowest level with the observables.
    for (size_t i = 0; i < obs_list.size(); i++) {
        std::set<sptr<raw_vertex_t>> obs;
        for (sptr<tanner::vertex_t> x : obs_list[i]) {
            obs.push_back(raw_connection_network.v_tanner_raw_map.at(x));
        }
        stim::simd_bits<SIMD_WIDTH> h(HASH_SIZE);
        h[i] = 1;

        curr_operator_tree_level.push_back(std::make_tuple(obs, h));
    }
    // Get relevant checks.
    std::vector<quop_t> checks;
    auto _checks = is_x_error ? tanner_graph.get_vertices_by_type(tanner::vertex_t::type::xparity)
                                : tanner_graph.get_vertices_by_type(tanner::vertex_t::type::zparity);
    for (size_t i = 0; i < _checks.size(); i++) {
        sptr<tanner::vertex_t> tpq = _checks[i];
        // Get the support of the check.
        std::set<sptr<raw_vertex_t>> stabilizer;
        for (sptr<tanner::vertex_t> tx : tanner_graph.get_neighbors(tpq)) {
            stabilizer.push_back(raw_connection_network.v_tanner_raw_map.at(tx));
        }
        stim::simd_bits<SIMD_WIDTH> h(HASH_SIZE);
        h[HASH_STABILIZER_OFFSET + i] = 1;

        checks.push_back(std::make_tuple(stabilizer, h));
    }
    // Perform a BFS via the operator tree.
    while (curr_operator_tree_level.size() && result.not_zero()) {
        // Check if any flags are a weight >2 error on any existing operators.
        std::vector<quop_t> next_level;
        for (quop_t op : curr_operator_tree_level) {
            auto op_qubits = std::get<0>(op);
            auto op_hash = std::get<1>(op);

            size_t k = 0;
            for (auto pair : flag_pairs) {
                sptr<raw_vertex_t> v1 = pair.first,
                                    v2 = pair.second;
                if (!result[k]) continue;   // This has already failed.
                if (op_qubits.count(v1) && op_qubits.count(v2)) {
                    // This is a weight-2 error.
                    result[k] = 0;
                }
            }
            // Now, expand the operator by multiplying it by each stabilizer.
            //
            // If the operator has already been multipled by a stabilizer (detectable via the
            // hash), we skip the product.
            //
            // If the operator and stabilizer share no qubits in common, then skip the product.
            for (quop_t stab : checks) {
                auto stab_qubits = std::get<0>(stab);
                auto stab_hash = std::get<1>(stab);

                stim::simd_bits<SIMD_WIDTH> common_hash_bits = op_hash & stab_hash;
                if (common_hash_bits.not_zero()) continue;

                std::set<raw_vertex_t> common_qubits = op_qubits * stab_qubits;
                if (common_qubits.empty()) continue;

                // Otherwise, everything is fine.
                std::set<raw_vertex_t> new_qubits = op_qubits + stab_qubits;
                stim::simd_bits<SIMD_WIDTH> new_hash = op_hash | stab_hash;
                next_level.push_back(std::make_tuple(new_qubits, new_hash));
            }
        }
        curr_operator_tree_level = std::move(next_level);
    }
    return result;
}

}   // protean
}   // qontra
