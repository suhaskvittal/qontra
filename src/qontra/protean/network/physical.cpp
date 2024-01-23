/* 
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */


#include "qontra/protean/network.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#include <PerfectMatching.h>

#include <algorithm>
#include <limits>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

using namespace vtils;

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

bool
PhysicalNetwork::join_qubits_with_identical_support() {
    bool mod = false;
    // These are vertices that are consumed by another.
    std::set<sptr<phys_vertex_t>> deleted_vertices;

    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<phys_vertex_t> pv = vertices[i];
        if (deleted_vertices.count(pv)) continue;
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
    
        auto pv_adj = get_neighbors(pv);
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
                mod = true;
            }
        }
    }
    // Delete the consumed vertices
    for (sptr<phys_vertex_t> x : deleted_vertices) delete_vertex(x);
    return mod;
}

bool
PhysicalNetwork::join_qubits_with_partial_support() {
    // Disable memoization.
    raw_connection_network.disable_memoization();

    bool mod = false;
    // These are vertices that are consumed by another.
    std::set<sptr<phys_vertex_t>> deleted_vertices;

    std::vector<sptr<phys_vertex_t>> _vertices(vertices);
    // Sort these vertices by degree.
    std::sort(_vertices.begin(), _vertices.end(),
            [&] (auto v, auto w) {
                return get_degree(v) > get_degree(w);
            });

    for (size_t i = 0; i < _vertices.size(); i++) {
        sptr<phys_vertex_t> pv = _vertices[i];
        if (deleted_vertices.count(pv)) continue;
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
    
        auto pv_adj = get_neighbors(pv);
        // Now, unlike identical_support, we will rank all qubits (that are non-data!!!)
        // by the amount of common neighbors. We will merge the best such qubit.
        sptr<phys_vertex_t> best_candidate = nullptr;
        size_t max_common_neighbors = 0;
        for (size_t j = i+1; j < _vertices.size(); j++) {
            sptr<phys_vertex_t> pw = _vertices[j];
            if (deleted_vertices.count(pw)) continue;
            if (pw->has_role_of_type(raw_vertex_t::type::data)) continue;

            size_t common_neighbors = 0;
            for (sptr<phys_vertex_t> x : pv_adj) {
                common_neighbors += contains(pw, x);
            }
            if (common_neighbors == get_degree(pw) && 
                    common_neighbors > max_common_neighbors) 
            {
                best_candidate = pw;
                max_common_neighbors = common_neighbors;
            }
        }
        if (best_candidate == nullptr) continue;
        consume(pv, best_candidate);
        deleted_vertices.insert(best_candidate);
        mod = true;
    }
    // Delete the consumed vertices
    for (sptr<phys_vertex_t> x : deleted_vertices) delete_vertex(x);
    return mod;
}

bool
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
                    insert_all(already_matched, {rv, rw});
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
        if (support.size() % 2 == 1) {
            already_matched.insert(support.back()); // It really hasn't, but whatever :)
            support.pop_back();
        }
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
            insert_all(visited, {rv, rw});
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
            role_to_phys[rfq] = pfq;
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
    return !all_proposed_flag_pairs.empty();
}

bool
PhysicalNetwork::add_connectivity_reducing_proxies() {
    bool mod = false;

    std::map<sptr<phys_vertex_t>, int> degree_update_map;
    std::vector<sptr<phys_edge_t>> new_edges;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        const size_t dg = get_degree(pv) + degree_update_map[pv];
        if (dg <= config.max_connectivity) continue;

        size_t slack_violation = dg - config.max_connectivity;
        std::vector<sptr<phys_vertex_t>> adj = get_neighbors(pv);
        // Also sort neighbors by degree.
        std::sort(adj.begin(), adj.end(),
                [&] (auto a, auto b)
                {
                    return get_degree(a)+degree_update_map[a] > get_degree(b)+degree_update_map[b];
                });

        sptr<phys_vertex_t> pprx = make_vertex();
        // Let k = slack_violation. Then the top k+1 neighbors will
        // be serviced by the proxy.
        std::vector<sptr<phys_vertex_t>> share_an_edge_with_pprx{pv};
        for (size_t i = 0; i < std::min(slack_violation+1, adj.size()); i++) {
            sptr<phys_vertex_t> px = adj[i];
            // Now, for each role in pv and px, add a proxy in the raw_connection_network.
            // Each added proxy will be a role for pprx.
            bool any_roles_added = false;
            for (sptr<raw_vertex_t> rv : pv->role_set) {
                for (sptr<raw_vertex_t> rx : px->role_set) {
                    sptr<raw_edge_t> re = raw_connection_network.get_edge(rv, rx);
                    if (re == nullptr) continue;
                    sptr<raw_vertex_t> rprx = raw_connection_network.add_proxy(re);

                    /*
                    std::set<uint64_t> good_ids{33, 231, 310, 130, 305, 222, 317, 70};
                    if (good_ids.count(rv->id) && good_ids.count(rx->id))
                    std::cout << "[ add_connectivity_reducing_proxies ] adding proxy " << print_v(rprx) 
                                << " between "
                                << print_v(pv) << "(" << print_v(rv) << ") and "
                                << print_v(px) << "(" << print_v(rx) << ")"
                                << " pointing [" << print_v(rv) << "] --> " << print_v(raw_connection_network.proxy_indirection_map[rprx][rv][0])
                                << " and [" << print_v(rx) << "] --> " << print_v(raw_connection_network.proxy_indirection_map[rprx][rx][0])
                                << std::endl;
                    */

                    // Get cycle of role (which is just the max of the cycles of rv and rx).
                    size_t cycle = std::max(pv->cycle_role_map.at(rv), px->cycle_role_map.at(rx));
                    pprx->push_back_role(rprx, cycle);
                    role_to_phys[rprx] = pprx;
                    any_roles_added = true;
                }
            }
            if (any_roles_added) {
                share_an_edge_with_pprx.push_back(px);
                // Delete px's edge with pv.
                sptr<phys_edge_t> pe = get_edge(pv, px);
                delete_edge(pe);
            }
        }
        if (pprx->role_set.empty()) {
            std::cerr << "[ add_connectivity_reducing_proxies ] found proxy " << print_v(pprx) << " with no roles\n";
            exit(1);
        }
        add_vertex(pprx);
        // Update the physical connectivity on our side.
        for (sptr<phys_vertex_t> x : share_an_edge_with_pprx) {
            new_edges.push_back(make_edge(pprx, x));
            degree_update_map[x]++;
        }
        mod = true;
    }
    for (auto e : new_edges) add_edge(e);
    return mod;
}

bool
PhysicalNetwork::contract_small_degree_qubits() {
    raw_connection_network.disable_memoization();

    bool mod = false;
    // Edges aren't added until after the loop ends. So, we need to track any updates.
    std::map<sptr<phys_vertex_t>, int> degree_update_map;
    std::vector<sptr<phys_edge_t>> new_edges;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
        
        const size_t dg = get_degree(pv) + degree_update_map[pv];
        if (dg > 2) continue;
        
        if (dg == 0) {
            // God forbid, I can't think of a case like this, so I'm printing out the
            // qubit.
            std::cerr << "[ warning ] found qubit with zero degree: " << print_v(pv) << "\n";
            std::cerr << "\troles =";
            for (auto r : pv->role_set) std::cerr << " " << print_v(r);
            std::cerr << "\n";
            delete_vertex(pv);
        } else if (dg == 1) {
            sptr<phys_vertex_t> px = get_neighbors(pv)[0];
            // If px is a data qubit, do not proceed.
            if (px->has_role_of_type(raw_vertex_t::type::data)) continue;

            // This is simple -- just have px consume pv and delete pv.
            consume(px, pv);
            delete_vertex(pv);
        } else {    // dg == 2
            sptr<phys_vertex_t> px1 = get_neighbors(pv)[0],
                                px2 = get_neighbors(pv)[1];
            // Simple: delete pv, and just add an edge between px1 and px2.
            // Move roles of pv to px1 or px2.
            if (!px1->has_role_of_type(raw_vertex_t::type::data)) {
                consume(px1, pv);
            } else if (!px2->has_role_of_type(raw_vertex_t::type::data)) {
                consume(px2, pv);
            } else {
                continue;
            }

            delete_vertex(pv);
            if (!contains(px1, px2)) {
                new_edges.push_back(make_edge(px1, px2));
                degree_update_map[px1]++;
                degree_update_map[px2]++;
            }
        }
        mod = true;
    }
    for (auto e : new_edges) add_edge(e);
    return mod;
}

bool
PhysicalNetwork::reallocate_edges() {
    // Top-down: trickle down the edges as much as possible.
    bool mod = false;
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
            mod |= test_and_move_edge_down(e);
        }
    }
    // Now, everything will have changed. Remove any layers with no edges.
    // A property of this function is that only the last layers should have no edges. So,
    // removing layers is rather straightforward:
    while (processor_layers.back().m() == 0 && processor_layers.size() > 1) {
        processor_layers.pop_back();
    }
    return mod;
}

stim::simd_bits<SIMD_WIDTH>
PhysicalNetwork::do_flags_protect_weight_two_error(
        std::set<std::pair<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>> flag_pairs,
        bool is_x_error) 
{
    // For the purposes of this algorithm, have result be active-low.
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
    const size_t OPERATOR_TREE_SIZE_LIMIT = 256L*1024L;
    const fp_t OPERATOR_WEIGHT_MULT = 2.0;

    // The first entry are the qubits in the operator.
    // THe second entry is a simple hash that allows us to check if two operators are equivalent.
    typedef std::tuple<std::set<sptr<raw_vertex_t>>, stim::simd_bits<SIMD_WIDTH>> quop_t;
    std::vector<quop_t> curr_operator_tree_level;

    const size_t HASH_SIZE = obs_list.size() + tanner_graph.get_checks().size()/2;
    const size_t HASH_STABILIZER_OFFSET = obs_list.size();
    // Populate the lowest level with the observables.
    size_t min_observable_size = std::numeric_limits<size_t>::max();
    for (size_t i = 0; i < obs_list.size(); i++) {
        std::set<sptr<raw_vertex_t>> obs;
        for (sptr<tanner::vertex_t> x : obs_list[i]) {
            obs.insert(raw_connection_network.v_tanner_raw_map.at(x));
        }
        min_observable_size = std::min(min_observable_size, obs.size());

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

        checks.emplace_back(stabilizer, h);
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
            if (op_qubits.size() == min_observable_size) {
                for (auto pair : flag_pairs) {
                    sptr<raw_vertex_t> v1 = pair.first,
                                        v2 = pair.second;
                    if (result[k] && op_qubits.count(v1) && op_qubits.count(v2)) {
                        // This is a weight-2 error.
                        result[k] = 0;
                    }
                    k++;
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
                std::set<sptr<raw_vertex_t>> new_qubits = op_qubits ^ stab_qubits;
                stim::simd_bits<SIMD_WIDTH> new_hash = op_hash | stab_hash;

                if (new_qubits.size() < OPERATOR_WEIGHT_MULT*min_observable_size) {
                    next_level.emplace_back(new_qubits, new_hash);
                }
            }
        }
        curr_operator_tree_level = std::move(next_level);
    }
    if (curr_operator_tree_level.size() > OPERATOR_TREE_SIZE_LIMIT) {
        std::cerr << "[ warning ] operator tree limit reached. exiting error propagation analysis.\n";
    }
    result.invert_bits();
    return result;
}

bool
PhysicalNetwork::recompute_cycle_role_maps() {
    raw_connection_network.disable_memoization();
    // Here, we want to enforce that for each check, each physical qubit has at most one role for that
    // check.
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        std::set<sptr<raw_vertex_t>> deleted_vertices;
        for (sptr<raw_vertex_t> r1 : pv->role_set) {
            for (sptr<raw_vertex_t> r2 : pv->role_set) {
                if (r1 <= r2) continue;
                // If r1 is deleted during a test_and_merge, break.
                if (deleted_vertices.count(r1)) break;
                if (deleted_vertices.count(r2)) continue;
                if (raw_connection_network.are_in_same_support(r1, r2) == nullptr) continue;

                auto deleted = raw_connection_network.merge(r1, r2);
                if (deleted != nullptr) deleted_vertices.insert(deleted);
            }
        }
        for (sptr<raw_vertex_t> r : deleted_vertices) {
            pv->delete_role(r);
            role_to_phys.erase(r);
        }
    }
    // Enable memoization.
    raw_connection_network.enable_memoization = true;
    // Basic idea: pack the roles together by parity check.
    //
    // We track an interaction graph, such that two checks share
    // an edge if any of their roles share a physical qubit. We
    // also track these conflicts.
    typedef std::tuple<sptr<raw_vertex_t>, sptr<raw_vertex_t>, sptr<phys_vertex_t>>
        conflict_t;

    struct int_v_t : public base::vertex_t {
        sptr<raw_vertex_t>              check;
        std::set<sptr<raw_vertex_t>>    support;
    };

    struct int_e_t : public base::edge_t {
        std::vector<conflict_t> conflicts;
    };

    TannerGraph& tanner_graph = raw_connection_network.tanner_graph;

    Graph<int_v_t, int_e_t> interaction_graph;
    size_t _id = 0;
    for (sptr<tanner::vertex_t> tpq : tanner_graph.get_checks()) {
        sptr<raw_vertex_t> rpq = raw_connection_network.v_tanner_raw_map.at(tpq);
        // Make interaction graph vertex.
        sptr<int_v_t> iv = interaction_graph.make_vertex();
        iv->id = _id++;
        iv->check = rpq;
        // Get all qubits involved in check.
        RawNetwork::parity_support_t& supp = raw_connection_network.get_support(rpq);
        iv->support = supp.all;

        interaction_graph.add_vertex(iv);
    }
    // Now, compute interaction graph edges.
    for (size_t i = 0; i < _id; i++) {
        sptr<int_v_t> iv = interaction_graph.get_vertex(i);
        for (size_t j = i+1; j < _id; j++) {
            sptr<int_v_t> iw = interaction_graph.get_vertex(j);

            sptr<int_e_t> ie = interaction_graph.make_edge(iv, iw);
            ie->is_undirected = true;
            for (sptr<raw_vertex_t> rx : iv->support) {
                if (rx->qubit_type == raw_vertex_t::type::data) continue;

                sptr<phys_vertex_t> px = role_to_phys[rx];
                for (sptr<raw_vertex_t> ry : iw->support) {
                    sptr<phys_vertex_t> py = role_to_phys[ry];
                    if (px == py) ie->conflicts.emplace_back(rx, ry, px);
                }
            }
            if (ie->conflicts.size()) interaction_graph.add_edge(ie);
        }
    }
    // Clear all the cycle role maps for each physical qubit.
    for (sptr<phys_vertex_t> pv : get_vertices()) pv->clear_roles();
    // Now, pack the roles. The idea here is for each iteration, schedule checks
    // that have no conflicts with each other.
    auto remaining_checks = interaction_graph.get_vertices();
    size_t cycle = 0;
    while (remaining_checks.size()) {
        std::vector<sptr<int_v_t>> checks_this_cycle;
        for (auto it = remaining_checks.begin(); it != remaining_checks.end(); ) {
            bool any_conflict = false;
            sptr<int_v_t> iv = *it;
            for (sptr<int_v_t> iw : checks_this_cycle) {
                if (interaction_graph.contains(iv, iw)) {
                    any_conflict = true;
                    break;
                }
            }
            if (any_conflict) {
                it++;
            } else {
                checks_this_cycle.push_back(iv);
                it = remaining_checks.erase(it);
            }
        }
        for (sptr<int_v_t> iv : checks_this_cycle) {
#ifdef PROTEAN_DEBUG
            std::cout << "support(" << print_v(iv->check) << "):";
#endif
            for (sptr<raw_vertex_t> rv : iv->support) {
                sptr<phys_vertex_t> pv = role_to_phys[rv];
                if (!pv->role_set.count(rv)) {
                    if (rv->qubit_type == raw_vertex_t::type::data) {
                        pv->push_back_role(rv, 0);
                    } else {
                        pv->push_back_role(rv, cycle);
                    }
                }
#ifdef PROTEAN_DEBUG
                std::cout << "\t" << print_v(rv) << " (" << print_v(pv)
                    << " @ " << pv->cycle_role_map.at(rv) << ")";
#endif
            }
#ifdef PROTEAN_DEBUG
            std::cout << "\n";
#endif
        }
        cycle++;
    }
    return true;
}

}   // protean
}   // qontra
