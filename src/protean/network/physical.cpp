/* author: Suhas Vittal date:   27 December 2023
 * */


#include "protean/network.h"

#include <defs/set_algebra.h>

#include <PerfectMatching.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

template <class V> using v_pair_t=std::pair<sptr<V>, sptr<V>>;

template <class T> inline std::pair<T, T>
make_ordered_pair(T x, T y) {
    if (x < y)  return std::make_pair(x, y);
    else        return std::make_pair(y, x);
}

PhysicalNetwork::PhysicalNetwork(TannerGraph& tgr)
    :Graph(),
    raw_connection_network(tgr),
    role_to_phys(),
    // Planarity tracking:
    processor_layers(),
    // Other tracking:
    occupied_tsvs(),
    // Other variables:
    id_ctr(0)
{
    push_back_new_processor_layer(); // The first ever processor layer.
    // Create a corresponding physical qubit for every vertex in raw_connection_network.
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        sptr<phys_vertex_t> pv = make_vertex();
        pv->push_back_role(rv);
        role_to_phys[rv] = pv;
        add_vertex(pv);
    }
    // And the corresponding edges as well.
    for (sptr<raw_edge_t> re : raw_connection_network.get_edges()) {
        sptr<raw_vertex_t> rsrc = raw_connection_network.get_source(re),
                            rdst = raw_connection_network.get_target(re);
        sptr<phys_edge_t> pe = make_edge(role_to_phys[rsrc], role_to_phys[rdst]);
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
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
    
        auto tmp = get_neighbors(pv);
        std::set<sptr<phys_vertex_t>> pv_adj(tmp.begin(), tmp.end());
        for (size_t j = i+1; j < vertices.size(); j++) {
            sptr<phys_vertex_t> pw = vertices[j];
            if (deleted_vertices.count(pw)) continue;
            if (pw->has_role_of_type(raw_vertex_t::type::data)) continue;

            size_t common_neighbors = 0;
            for (sptr<phys_vertex_t> x : pv_adj) {
                common_neighbors += contains(pw, x);
            }
            if (common_neighbors == pv_adj.size()) {
                consume(pv, pw); // Now, pv contains all the roles of pw.
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

    // The entries should be ordered pairs to avoid double counting.
    typedef std::set<v_pair_t<raw_vertex_t>>    flag_pair_set_t;

    std::map<sptr<raw_vertex_t>, flag_pair_set_t> proposed_flag_pair_map;
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

        // Perform the perfect matching:
        typedef std::tuple<sptr<raw_vertex_t>, sptr<raw_vertex_t>, size_t>  weighted_edge_t;
        std::vector<weighted_edge_t> edge_list;
        size_t max_edge_weight = 0;
        // In case anything is already matched, we will just forego including them in the matching
        // problem and pre-assign them.
        flag_pair_set_t flag_pairs;

        std::set<sptr<raw_vertex_t>> visited;
        std::set<sptr<raw_vertex_t>> already_matched;
        for (auto ita = support.begin(); ita != support.end(); ) {
edge_build_outer_loop_start:
            if (ita == support.end()) break;
            sptr<raw_vertex_t> rv = *ita;
            // We will the physical qubit (specifically its neighbors) each tick of the
            // inner loop.
            // Note that we only care about neighbors that are parity qubits.
            sptr<phys_vertex_t> pv = role_to_phys[rv];
            std::set<sptr<phys_vertex_t>> p_neighbors;
            for (sptr<phys_vertex_t> x : get_neighbors(pv)) {
                if (x->has_role_of_type(raw_vertex_t::type::xparity)
                    || x->has_role_of_type(raw_vertex_t::type::zparity))
                {
                    p_neighbors.insert(x);
                }
            }
            for (auto itb = ita+1; itb != support.end(); ) {
                sptr<raw_vertex_t> rw = *itb;
            
                if (all_proposed_flag_pairs.count(make_ordered_pair(rv, rw))) {
                    // Then, rv and rw are already assigned. Move on.
                    // Erase the vertices from the support.
                    // 
                    // Unsafe implementation hack: we need to delete two iterators,
                    // but erase does not guarantee that ita/itb will be valid after
                    // the erase of one. But, logically deleting an element later
                    // in the array (itb) should not affect an earlier element (ita).
                    itb = support.erase(itb);
                    ita = support.erase(ita);
                    // Add this pair to the flag pairs.
                    flag_pairs.insert(make_ordered_pair(rv, rw));
                    already_matched.insert(rv);
                    already_matched.insert(rw);
                    goto edge_build_outer_loop_start;
                }
                // Get corresponding physical vertices, get number of common neighbors,
                // and set the weight of edge to that.
                sptr<phys_vertex_t> pw = role_to_phys[rw];
                size_t common_neighbors = 0;
                for (sptr<phys_vertex_t> x : p_neighbors) {
                    common_neighbors += contains(pw, x);
                }
                edge_list.emplace_back(rv, rw, common_neighbors);
                max_edge_weight = std::max(max_edge_weight, common_neighbors);
                itb++;
            }
            ita++;
        }
        // Now, create the MWPM object.
        PerfectMatching pm(support.size(), edge_list.size());
        pm.options.verbose = false;
        // Assign each node in the support to some node id.
        BijectiveMap<sptr<raw_vertex_t>, size_t> data_node_map;
        for (size_t i = 0; i < support.size(); i++) data_node_map.put(i, support[i]);
        // Get the edges from the edge_list. Note that some of these edges are bad (i.e.
        // containing a pre-matched vertex).
        //
        // We also need to reweight the edge by max_edge_weight - W -- this is a MWPM solver.
        for (weighted_edge_t we : edge_list) {
            sptr<raw_vertex_t> v = std::get<0>(we),
                               w = std::get<1>(we);
            size_t weight = max_edge_weight - std::get<2>(we) + 1;
            if (already_matched.count(v) || already_matched.count(w)) continue;
            pm.AddEdge(data_node_map.at(v), data_node_map.at(w), weight);
        }
        pm.Solve();
        // Get the PM result.
        for (sptr<raw_vertex_t> rv : support) {
            if (visited.count(rv)) continue;

            size_t i = data_node_map.at(rv);
            size_t j = pm.GetMatch(i);
            sptr<raw_vertex_t> rw = data_node_map.at(j);
            
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

    std::cout << "[ debug ] flags before error-propagation analysis: " << all_proposed_flag_pairs.size() << "\n";
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
    std::cout << "[ debug ] flags after error-propagation analysis: " << all_proposed_flag_pairs.size() << "\n";

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

        sptr<phys_vertex_t> pfq = make_vertex();
        std::set<sptr<phys_vertex_t>> qubits_connected_to_pfq{pv1, pv2};
        for (auto& p2 : role_set) {
            sptr<raw_vertex_t> rpq = p2.first,
                                rfq = p2.second;
            sptr<phys_vertex_t> ppq = role_to_phys[rpq];

            // Compute cycle for role. This is the maximum of the cycles of rv1, rv2, and rfq.
            size_t cycle = std::max({pv1->cycle_role_map.at(rv1),
                                    pv2->cycle_role_map.at(rv2),
                                    ppq->cycle_role_map.at(rpq)});
            pfq->push_back_role(rfq, cycle);
            qubits_connected_to_pfq.insert(ppq);
            // Delete edge between pv1/pv2 and ppq (if it exists).
            if (contains(pv1, ppq)) {
                delete_edge(get_edge(pv1, ppq));
            } 
            if (contains(pv2, ppq)) {
                delete_edge(get_edge(pv2, ppq));
            }
        }
        // Add edges for pfq.
        add_vertex(pfq);
        for (sptr<phys_vertex_t> x : qubits_connected_to_pfq) {
            sptr<phys_edge_t> pe = make_edge(pfq, x);
            add_edge(pe);
        }
    }
    reallocate_edges();
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

        sptr<phys_vertex_t> pprx = make_vertex();
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
                    pprx->push_back_role(rprx, cycle);
                }
            }
            share_an_edge_with_pprx.insert(px);
            // Delete px's edge with pv.
            sptr<phys_edge_t> pe = get_edge(pv, px);
            delete_edge(pe);
        }
        // Update the physical connectivity on our side.
        for (sptr<phys_vertex_t> x : share_an_edge_with_pprx) {
            sptr<phys_edge_t> e = make_edge(pprx, x);
            add_edge(e);
        }
    }
    // If any proxies were added, we need to check the connectivity of these proxies.
    if (were_any_proxies_added) add_connectivity_reducing_proxies();
}

void
PhysicalNetwork::contract_small_degree_qubits() {
    std::vector<sptr<phys_edge_t>> new_edges;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
        
        const size_t dg = get_degree(pv);
        if (dg > 2) continue;
        
        if (dg == 1) {
            sptr<phys_vertex_t> px = get_neighbors(pv)[0];
            if (px->has_role_of_type(raw_vertex_t::type::data)) continue;

            std::cout << "[ debug ] contracting " << print_v(pv) << " into " << print_v(px) << "\n";

            // This is simple -- just have px consume pv and delete pv.
            consume(px, pv);
            delete_vertex(pv);
        } else {    // dg == 2
            sptr<phys_vertex_t> px1 = get_neighbors(pv)[0],
                                px2 = get_neighbors(pv)[1];
            // Simple: delete pv, and just add an edge between px1 and px2.
            //
            // Move roles of pv to px1 or px2.
            if (!px1->has_role_of_type(raw_vertex_t::type::data)) {
                consume(px1, pv);
            } else if (!px2->has_role_of_type(raw_vertex_t::type::data)) {
                consume(px2, pv); 
            } else {
                continue;
            }
            std::cout << "[ debug ] contracting " << print_v(pv) << " into u(" 
                    << print_v(px1) << "," << print_v(px2) << ")\n";

            delete_vertex(pv);
            if (!contains(px1, px2)) {
                new_edges.push_back(make_edge(px1, px2));
            }
        }
    }
    for (auto e : new_edges) {
        add_edge(e);
        std::cout << "[ debug ] adding edge " << print_e<phys_vertex_t>(e) << " to L" << e->tsv_layer << "\n";
    }
    reallocate_edges();
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
    typedef std::tuple<std::set<sptr<raw_vertex_t>>, stim::simd_bits<SIMD_WIDTH>> quop_t;
    std::vector<quop_t> curr_operator_tree_level;

    const size_t HASH_SIZE = obs_list.size() + tanner_graph.get_checks().size()/2;
    const size_t HASH_STABILIZER_OFFSET = obs_list.size();
    // Populate the lowest level with the observables.
    for (size_t i = 0; i < obs_list.size(); i++) {
        std::set<sptr<raw_vertex_t>> obs;
        for (sptr<tanner::vertex_t> x : obs_list[i]) {
            obs.insert(raw_connection_network.v_tanner_raw_map.at(x));
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
            stabilizer.insert(raw_connection_network.v_tanner_raw_map.at(tx));
        }
        stim::simd_bits<SIMD_WIDTH> h(HASH_SIZE);
        h[HASH_STABILIZER_OFFSET + i] = 1;

        checks.push_back(std::make_tuple(stabilizer, h));
    }
    // Perform a BFS via the operator tree.
    while (curr_operator_tree_level.size() > 0
        && curr_operator_tree_level.size() < OPERATOR_TREE_SIZE_LIMIT 
        && result.not_zero()) 
    {
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

                std::set<sptr<raw_vertex_t>> common_qubits = op_qubits * stab_qubits;
                if (common_qubits.empty()) continue;

                // Otherwise, everything is fine.
                std::set<sptr<raw_vertex_t>> new_qubits = op_qubits + stab_qubits;
                stim::simd_bits<SIMD_WIDTH> new_hash = op_hash | stab_hash;
                next_level.push_back(std::make_tuple(new_qubits, new_hash));
            }
        }
        curr_operator_tree_level = std::move(next_level);
    }
    if (curr_operator_tree_level.size() > OPERATOR_TREE_SIZE_LIMIT) {
        std::cerr << "[ status ] operator tree limit reached. exiting error propagation analysis.\n";
    }
    return result;
}

schedule_t
PhysicalNetwork::make_schedule() {
    return schedule_t();
}

void
PhysicalNetwork::reallocate_edges() {
    std::cout << "[ debug ] reallocating edges | # of layers = " << get_thickness() << "\n";
    // Top-down: trickle down the edges as much as possible.
    bool any_edges_were_reallocated = false;
    for (size_t k = get_thickness()-1; k > 0; k--) {
        ProcessorLayer& curr = processor_layers[k];
        ProcessorLayer& next = processor_layers[k-1];
        // Get a sorted list of edges in curr, sorted from least to greatest
        // by their max_endpoint_degree in next.
        auto edge_list = curr.get_edges();
        std::sort(edge_list.begin(), edge_list.end(),
                [&] (auto e, auto f)
                {
                    return next.get_max_endpoint_degree(e) < next.get_max_endpoint_degree(f);
                });
        // Now, try to put as many edges into next.
        for (sptr<phys_edge_t> e : edge_list) {
            any_edges_were_reallocated |= test_and_move_edge_down(e);
        }
    }
    // Now, everything will have changed. Remove any layers with no edges.
    // A property of this function is that only the last layers should have no edges. So,
    // removing layers is rather straightforward:
    while (processor_layers.back().m() == 0 && processor_layers.size() > 1) {
        processor_layers.pop_back();
    }
    if (any_edges_were_reallocated) reallocate_edges();
}

}   // protean
}   // qontra
