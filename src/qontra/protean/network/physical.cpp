/* 
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "qontra/protean/network.h"
#include "qontra/protean/scheduler.h"
#include "qontra/protean/yield_sim.h"

#include <qontra/graph/algorithms/coloring.h>

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#ifdef PROTEAN_PERF
#include <vtils/timer.h>
#endif

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

PhysicalNetwork::PhysicalNetwork(TannerGraph* tgr)
    :Graph(),
    tanner_graph(tgr),
    raw_connection_network(std::make_unique<RawNetwork>(tgr)),
    role_to_phys(),
    // Planarity tracking:
    processor_layers(),
    // Other tracking:
    check_color_map(),
    // Other variables:
    id_ctr(0),
    round_latency(0),
    round_cnots(0),
    expected_collisions(0),
    x_memory(),
    z_memory()
{
    push_back_new_processor_layer(); // The first ever processor layer.
    // Create a corresponding physical qubit for every vertex in raw_connection_network.
    for (sptr<raw_vertex_t> rv : raw_connection_network->get_vertices()) {
        sptr<phys_vertex_t> pv = make_and_add_vertex();
        pv->push_back_role(rv);
        role_to_phys[rv] = pv;
    }
    // And the corresponding edges as well.
    for (sptr<raw_edge_t> re : raw_connection_network->get_edges()) {
        sptr<raw_vertex_t> rsrc = re->get_source<raw_vertex_t>(),
                            rdst = re->get_target<raw_vertex_t>();
        make_and_add_edge(role_to_phys.at(rsrc), role_to_phys.at(rdst));
    }
    raw_connection_network->disable_memoization();
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
    typedef std::set<v_pair_t<raw_vertex_t>>    flag_pair_set_t;

    // Track flag pairs by check role.
    std::map<sptr<raw_vertex_t>, flag_pair_set_t> proposed_flag_pair_map;
    // Because we want minimize resource usage, also track the flags incident upon a physical qubit.
    // Checks on the same qubit will use the same flag pairs.
    std::map<sptr<phys_vertex_t>, flag_pair_set_t> phys_proposed_flag_pair_map;
    // Finally, track all flag pairs.
    flag_pair_set_t all_proposed_flag_pairs;
    flag_pair_set_t x_flags, z_flags;

#ifdef PROTEAN_PERF
    Timer timer;
    fp_t t;
    timer.clk_start();
#endif

    for (sptr<raw_vertex_t> rpq : raw_connection_network->get_vertices()) {
        if (!rpq->is_check()) continue;
        // So, first we will use _support, which is the set of data qubits. Then, we will
        // convert it to a vector.
        std::set<sptr<raw_vertex_t>> _support(raw_connection_network->get_support(rpq).data);

        flag_pair_set_t flag_pairs;
        std::set<sptr<raw_vertex_t>> already_matched;
        // Check if rpq's physical qubit has another role that has already given flags in
        // the proposed_flag_pair_map.
        sptr<phys_vertex_t> ppq = role_to_phys.at(rpq);
        for (const auto& fp : phys_proposed_flag_pair_map[ppq]) {
            sptr<raw_vertex_t> rv = fp.first,
                                rw = fp.second;
            if (_support.count(rv) && _support.count(rw)) {
                flag_pairs.insert(fp);
                insert_all(already_matched, {rv, rw});
                _support.erase(rv);
                _support.erase(rw);
            }
        }

        std::vector<sptr<raw_vertex_t>> support(_support.begin(), _support.end());

        // Perform the perfect matching:
        typedef std::tuple<sptr<raw_vertex_t>, sptr<raw_vertex_t>, size_t>  weighted_edge_t;
        std::vector<weighted_edge_t> edge_list;
        size_t max_edge_weight = 0;
        // In case anything is already matched, we will just forego including them in the matching
        // problem and pre-assign them.
        for (auto ita = support.begin(); ita != support.end(); ) {
edge_build_outer_loop_start:
            if (ita == support.end()) break;
            sptr<raw_vertex_t> rv = *ita;
            // We will the physical qubit (specifically its neighbors) each tick of the
            // inner loop.
            // Note that we only care about neighbors that are parity qubits.
            std::set<sptr<phys_vertex_t>> p_neighbors;
            for (sptr<raw_vertex_t> rx : raw_connection_network->get_neighbors(rv)) {
                if (rx->is_check()) {
                    p_neighbors.insert(role_to_phys.at(rx));
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
                sptr<phys_vertex_t> pw = role_to_phys.at(rw);
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
        for (size_t i = 0; i < support.size(); i++) {
            data_node_map.put(i, support[i]);
        }
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
            size_t i = data_node_map.at(rv);
            size_t j = pm.GetMatch(i);
            if (i > j) continue;
            sptr<raw_vertex_t> rw = data_node_map.at(j);
            
            flag_pairs.insert(make_ordered_pair(rv, rw));
        }
        // Update the flag proposal structures.
        insert_range(phys_proposed_flag_pair_map[ppq], flag_pairs);
        insert_range(all_proposed_flag_pairs, flag_pairs);
        proposed_flag_pair_map[rpq] = flag_pairs;
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            insert_range(x_flags, flag_pairs);
        } else {
            insert_range(z_flags, flag_pairs);
        }
    }
#ifdef PROTEAN_PERF
    t = timer.clk_end();
    std::cout << "[ make_flags ] flag identification took " << t*1e-9 << "s" << std::endl;
    timer.clk_start();
#endif
    // Now, we must run an error-propagation analysis to identify any useless flags.
    std::array<flag_pair_set_t, 2> flag_pair_sets{x_flags, z_flags};
    if (config.enable_flag_reduction) {
        size_t flags_removed = 0;
        const size_t total_flags = x_flags.size() + z_flags.size();
        for (size_t i = 0; i < 2; i++) {
            flag_pair_set_t& flags = flag_pair_sets[i];
            stim::simd_bits<SIMD_WIDTH> indicator_bits = do_flags_protect_weight_two_error(flags, i==0);
            // Remove any flags whose indicator bit is 0.
            size_t k = 0;
            for (auto it = flags.begin(); it != flags.end(); ) {
                if (indicator_bits[k]) {
                    it++;
                } else {
                    all_proposed_flag_pairs -= *it;
                    it = flags.erase(it);
                    flags_removed++;
                }
                k++;
            }
        }
        std::cout << "[ status ] removed flags = " << flags_removed << " of " << total_flags << std::endl;
#ifdef PROTEAN_PERF
        t = timer.clk_end();
        std::cout << "[ make_flags ] flag removal analysis took " << t*1e-9 << "s" << std::endl;
        timer.clk_start();
#endif
    }
    // First, we will create flags for each parity qubit. This ignores the case where (q1, q2)
    // is a flag in two different parity qubits. However, note we are only interested in establishing
    // roles here. These common roles will eventually be compacted into one flag.
    //
    // We track common roles in the following structure. The values are a parity-qubit and flag-qubit pair.
    std::map<v_pair_t<raw_vertex_t>, std::set<v_pair_t<raw_vertex_t>>>
        flag_pairs_to_flag_roles;

    for (auto& [ rpq, flags ] : proposed_flag_pair_map) {
        for (auto it = flags.begin(); it != flags.end(); ) {
            auto& fp = *it;
            if (all_proposed_flag_pairs.count(fp)) {
                sptr<raw_vertex_t> rv = fp.first,
                                    rw = fp.second;
                sptr<raw_vertex_t> rfq = raw_connection_network->add_flag(rv, rw, rpq);
                // Update flag_pair_to_flag_roles
                flag_pairs_to_flag_roles[fp].insert(std::make_pair(rpq, rfq));
                it++;
            } else {
                it = flags.erase(it);
            }
        }
    }
    // Now, we can make the physical qubits and update the connectivity.
    std::set<sptr<phys_edge_t>> deleted_edges;
    std::set<sptr<phys_edge_t>> new_edges;
    for (auto& [ rv_pair, role_set ] : flag_pairs_to_flag_roles) {
        auto [ rv, rw ] = rv_pair;
        sptr<phys_vertex_t> pv = role_to_phys.at(rv),
                            pw = role_to_phys.at(rw);
        // Make a separate physical qubit for the X and Z flags.
        for (int i = 0; i < (config.force_xz_flag_merge ? 1 : 2); i++) {
            sptr<phys_vertex_t> pfq;;
            if (!config.force_unopt_flags) {
                pfq = make_vertex(id_ctr++);
            }
            // Track connected qubits and deleted edges.
            std::set<sptr<phys_edge_t>> local_deleted_edges;
            std::set<sptr<phys_vertex_t>> qubits_connected_to_pfq{pv, pw};
            for (auto& [ rpq, rfq ] : role_set) {
                if (!config.force_xz_flag_merge && 
                        ((i == 0) ^ (rpq->qubit_type == raw_vertex_t::type::xparity))) 
                {
                    continue;
                }
                if (config.force_unopt_flags) {
                    pfq = make_and_add_vertex();
                }
                sptr<phys_vertex_t> ppq = role_to_phys.at(rpq);
                // Compute cycle for role. This is the maximum of the cycles of rv1, rv2, and rfq.
                size_t cycle = std::max({pv->cycle_role_map.at(rv),
                                        pw->cycle_role_map.at(rw),
                                        ppq->cycle_role_map.at(rpq)});
                pfq->push_back_role(rfq, cycle);
                role_to_phys[rfq] = pfq;
                if (config.force_unopt_flags) {
                    // Add the edges now and immediately update deleted_edges.
                    insert_all(new_edges, { make_edge(pfq, pv), make_edge(pfq, pw), make_edge(pfq, ppq) });
                    insert_all(deleted_edges, { get_edge(pv, ppq), get_edge(pw, ppq) });
                } else {
                    qubits_connected_to_pfq.insert(ppq);
                    // Delete edge between pv1/pv2 and ppq (if it exists).
                    insert_all(local_deleted_edges, { get_edge(pv, ppq), get_edge(pw, ppq) });
                }
            }
            // Add edges for pfq.
            if (!config.force_unopt_flags && qubits_connected_to_pfq.size() > 2) {
                add_vertex(pfq);
                insert_range(deleted_edges, local_deleted_edges);
                for (sptr<phys_vertex_t> x : qubits_connected_to_pfq) {
                    new_edges.insert(make_edge(pfq, x));
                }
            }
        }
    }
    for (sptr<phys_edge_t> e : deleted_edges) delete_edge(e);
    for (sptr<phys_edge_t> e : new_edges) add_edge(e);
#ifdef PROTEAN_PERF
    t = timer.clk_end();
    std::cout << "[ make_flags ] flag creation and placement took " << t*1e-9 << "s" << std::endl;
#endif
    return !all_proposed_flag_pairs.empty();
}

bool
PhysicalNetwork::add_connectivity_reducing_proxies() {
    std::cout << "[ status ] qubits = " << n() 
        << ", max connectivity = " << get_max_degree()
        << ", mean connectivity = " << get_mean_degree() << std::endl;
    bool mod = false;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        const size_t dg = get_degree(pv);
        if (dg <= config.max_connectivity) continue;

        int slack_violation = dg - config.max_connectivity;
        std::vector<sptr<phys_vertex_t>> adj = get_neighbors(pv);
        // Also sort neighbors by degree.
        std::sort(adj.begin(), adj.end(),
                [&] (auto a, auto b)
                {
                    return get_degree(a) > get_degree(b);
                });
        sptr<phys_vertex_t> pprx = make_and_add_vertex();
        // Let k = slack_violation. Then the top k+1 neighbors will
        // be serviced by the proxy.
        std::set<sptr<phys_vertex_t>> visited;
        std::vector<sptr<phys_vertex_t>> share_an_edge_with_pprx{pv};
        std::vector<sptr<phys_edge_t>> deleted_edges;
        slack_violation++;
        // We want to find triangles amongst the neighbors (i.e. a triple (pv, px, py)). This
        // means that the proxies reduce the connectivity of three vertices.
        bool searching_for_triangle = false;
        sptr<phys_vertex_t> pw = nullptr;

        size_t i = 0;
        while (slack_violation > 0 && visited.size() < adj.size()) {
            if (i == adj.size()) i = 0;
            sptr<phys_vertex_t> px = adj[i];
            if (config.enable_proxy_triangle_search && searching_for_triangle && px == pw) {
                // We could not find a triangle.
                i++;
                searching_for_triangle = false;
                continue;
            }
            if (visited.count(px)) {
                i++;
                continue;
            }
            if (config.enable_proxy_triangle_search && searching_for_triangle && !contains(pw, px)) {
                i++;
                continue;
            }
            // Now, for each role in pv and px, add a proxy in the raw_connection_network.
            // Each added proxy will be a role for pprx.
            bool any_roles_added = false;
            for (sptr<raw_vertex_t> rv : pv->role_set) {
                for (sptr<raw_vertex_t> rx : px->role_set) {
                    sptr<raw_edge_t> re = raw_connection_network->get_edge(rv, rx);
                    if (re == nullptr) continue;
                    sptr<raw_vertex_t> rprx = raw_connection_network->add_proxy(re);

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
                deleted_edges.push_back(get_edge(pv, px));
                slack_violation--;
                // Update searching_for_triangle
                if (config.enable_proxy_triangle_search) {
                    if (searching_for_triangle) {
                        deleted_edges.push_back(get_edge(pw, px));
                        searching_for_triangle = false;
                        std::cout << "found triangle" << std::endl;
                    } else {
                        searching_for_triangle = true;
                        pw = px;
                    }
                }
            }
            visited.insert(px);
            i++;
        }
        // Update the physical connectivity on our side.
        if (share_an_edge_with_pprx.size() > 2) {
            for (sptr<phys_edge_t> _pe : deleted_edges) {
                delete_edge(_pe);
            }
            for (sptr<phys_vertex_t> x : share_an_edge_with_pprx) {
                make_and_add_edge(pprx, x);
            }
            mod = true;
        } else {
            std::cerr << "TODO: handle proxy delete!" << std::endl;
            exit(1);
        }
    }
    return mod;
}

bool
PhysicalNetwork::contract_small_degree_qubits() {
    bool mod = false;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        if (pv->has_role_of_type(raw_vertex_t::type::data)) continue;
        
        const size_t dg = get_degree(pv);
        if (dg >= config.max_connectivity) continue;
        
        if (dg == 0) {
            // God forbid, I can't think of a case like this, so I'm exiting.
            std::cerr << "[ warning ] found qubit with zero degree: " << print_v(pv) << "\n";
            std::cerr << "\troles =";
            for (auto r : pv->role_set) std::cerr << " " << print_v(r);
            std::cerr << std::endl;
            exit(1);
        } else if (dg == 1) {
            sptr<phys_vertex_t> px = get_neighbors(pv)[0];
            // If px is a data qubit, do not proceed.
            if (px->has_role_of_type(raw_vertex_t::type::data)) continue;

            // This is simple -- just have px consume pv and delete pv.
            consume(px, pv);
            delete_vertex(pv);
        } else if (dg == 2) {
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
                make_and_add_edge(px1, px2);
            }
        } else {
            // Otherwise, this is some vertex with degree less than the max connectivity. If any
            // of its neighbors have sufficient slack to absorb it, we will do so.
            auto adj = get_neighbors(pv);
            for (sptr<phys_vertex_t> px : adj) {
                int slack = dg + get_degree(px) - 1;
                if (px->has_role_of_type(raw_vertex_t::type::data)
                        || slack > config.max_connectivity) continue;
                consume(px, pv);
                delete_vertex(pv);
                // Update neighbors of pv.
                for (sptr<phys_vertex_t> py : adj) {
                    if (px == py) continue;
                    if (!contains(px, py)) {
                        make_and_add_edge(px, py);
                    }
                }
                goto did_mod_network;
            }
            continue;
        }
did_mod_network:
        mod = true;
    }
    return mod;
}

bool
PhysicalNetwork::reallocate_edges() {
    if (!graph_has_changed) return false;
    // Top-down: trickle down the edges as much as possible.
    bool mod = false;
    for (size_t k = get_thickness()-1; k > 0; k--) {
        uptr<ProcessorLayer>& curr = processor_layers[k],
                            & next = processor_layers[k-1];
        // Get a sorted list of edges in curr, sorted from least to greatest
        // by their max_endpoint_degree in next.
        auto edge_list = curr->get_edges();
        std::sort(edge_list.begin(), edge_list.end(),
                [&] (auto e, auto f)
                {
                    return next->get_max_endpoint_degree(e) < next->get_max_endpoint_degree(f);
                });
        // Now, try to put as many edges into next.
        for (sptr<phys_edge_t> e : edge_list) {
            mod |= test_and_move_edge_down(e);
        }
    }
    // Now, everything will have changed. Remove any layers with no edges.
    // A property of this function is that only the last layers should have no edges. So,
    // removing layers is rather straightforward:
    while (processor_layers.back()->m() == 0 && processor_layers.size() > 1) {
        processor_layers.pop_back();
    }
    return mod;
}

stim::simd_bits<SIMD_WIDTH>
PhysicalNetwork::do_flags_protect_weight_two_error(
        std::set<std::pair<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>> flag_pairs,
        bool is_x_error) 
{
    // For the purposes of this algorithm, have the bitvectors be active-low.
    stim::simd_bits<SIMD_WIDTH> is_weight_two(flag_pairs.size()),   // 0 = good
                                is_stabilizer(flag_pairs.size());   // 1 = good
    is_weight_two.invert_bits();
    is_stabilizer.invert_bits();
    // Algorithm: we need to check against any logical operator. This can be rather exhaustive.
    // We will use the concept of an operator tree.
    // 
    // The lowest level of the tree corresponds to the observables of the code.
    auto obs_list = tanner_graph->get_obs(is_x_error);

    // We will only examine a single level of the tree at any time to avoid high memory overheads of
    // loading the entire tree.
    const size_t OPERATOR_TREE_SIZE_LIMIT = 4L*1024L;

    // The first entry are the qubits in the operator.
    // THe second entry is a simple hash that allows us to check if two operators are equivalent.
    typedef std::tuple<std::set<sptr<raw_vertex_t>>, stim::simd_bits<SIMD_WIDTH>> quop_t;
    std::vector<quop_t> curr_operator_tree_level;

    const size_t HASH_SIZE = obs_list.size() + tanner_graph->get_checks().size()/2;
    const size_t HASH_STABILIZER_OFFSET = obs_list.size();
    // Populate the lowest level with the observables.
    for (size_t i = 0; i < obs_list.size(); i++) {
        std::set<sptr<raw_vertex_t>> obs;
        for (sptr<tanner::vertex_t> x : obs_list[i]) {
            obs.insert(raw_connection_network->v_tanner_raw_map.at(x));
        }
        stim::simd_bits<SIMD_WIDTH> h(HASH_SIZE);
        h[i] = 1;
        curr_operator_tree_level.push_back(std::make_tuple(obs, h));
    }
    // Get relevant checks.
    std::vector<quop_t> checks;
    auto _checks = is_x_error ? tanner_graph->get_vertices_by_type(tanner::vertex_t::type::xparity)
                                : tanner_graph->get_vertices_by_type(tanner::vertex_t::type::zparity);
    for (size_t i = 0; i < _checks.size(); i++) {
        sptr<tanner::vertex_t> tpq = _checks[i];
        // Get the support of the check.
        std::set<sptr<raw_vertex_t>> stabilizer;
        for (sptr<tanner::vertex_t> tx : tanner_graph->get_neighbors(tpq)) {
            stabilizer.insert(raw_connection_network->v_tanner_raw_map.at(tx));
        }
        stim::simd_bits<SIMD_WIDTH> h(HASH_SIZE);
        h[HASH_STABILIZER_OFFSET + i] = 1;

        checks.emplace_back(stabilizer, h);
    }
    // Perform a BFS via the operator tree.
    stim::simd_bits<SIMD_WIDTH> result(flag_pairs.size());
    result.invert_bits();
    do {
        // Check if any flags are a weight >2 error on any existing operators.
        std::vector<quop_t> next_level;
        for (quop_t op : curr_operator_tree_level) {
            auto op_qubits = std::get<0>(op);
            auto op_hash = std::get<1>(op);

            size_t k = 0;
            for (const auto& [v1, v2] : flag_pairs) {
                if (result[k]) {
                    if (op_qubits.count(v1) && op_qubits.count(v2) && op_qubits.size() <= 12) {
                        // This is a weight-2 error.
                        is_weight_two[k] = 0;
                    }
                }
                // Check if op | v1 | v2 forms a stabilizer.
                std::set<sptr<raw_vertex_t>> op_v1_v2 = op_qubits ^ std::set<sptr<raw_vertex_t>>{v1, v2};
                for (quop_t ch : checks) {
                    if (std::get<0>(ch) == op_v1_v2) is_stabilizer[k] = 0;
                }
                k++;
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
                if (new_qubits.size() < 12) {
                    next_level.emplace_back(new_qubits, new_hash);
                }
            }
        }
        curr_operator_tree_level = std::move(next_level);

        result = is_weight_two & is_stabilizer;
    } while (curr_operator_tree_level.size() > 0
        && curr_operator_tree_level.size() <= OPERATOR_TREE_SIZE_LIMIT 
        && result.not_zero());
    if (curr_operator_tree_level.size() > OPERATOR_TREE_SIZE_LIMIT) {
        std::cerr << "[ warning ] operator tree limit reached. exiting error propagation analysis.\n";
    }
    // We want all flags where is_weight_two = 0 and is_stabilizer = 1.
    is_weight_two.invert_bits();
    result = is_weight_two & is_stabilizer;
    return result;
}

bool
PhysicalNetwork::recompute_cycle_role_maps() {
    // Here, we want to enforce that for each check, each physical qubit has at most one role for that
    // check.
    //
    // We also want to compress the flag roles, such that any physical qubit has
    // at most 1 Z-flag and at most 1 X-flag.
#ifdef PROTEAN_PERF
    Timer timer;
    timer.clk_start();
#endif
    for (sptr<raw_vertex_t> rpq : raw_connection_network->get_vertices()) {
        if (!rpq->is_check()) continue;
        auto& support = raw_connection_network->get_support(rpq);

        std::set<sptr<raw_vertex_t>> deleted_vertices;
        for (sptr<raw_vertex_t> rx : support.all) {
            if (deleted_vertices.count(rx)) continue;
            sptr<phys_vertex_t> px = role_to_phys.at(rx);
            for (sptr<raw_vertex_t> ry : px->role_set) {
                if (rx == ry) continue;
                if (deleted_vertices.count(ry)) continue;

                bool rx_ry_are_same_flag =
                    (rx->qubit_type == raw_vertex_t::type::flag
                        && ry->qubit_type == raw_vertex_t::type::flag)
                    &&
                    (raw_connection_network->x_flag_set.count(rx)
                        == raw_connection_network->x_flag_set.count(ry));
                // Try and merge rx and ry.
                if (support.all.count(ry) || rx_ry_are_same_flag) {
                    auto deleted = raw_connection_network->merge(rx, ry);
                    if (deleted != nullptr) {
                        deleted_vertices.insert(deleted);
                        if (deleted == rx) break;
                    }
                }
            }
        }
        for (sptr<raw_vertex_t> r : deleted_vertices) {
            sptr<phys_vertex_t> p = role_to_phys.at(r);
            p->delete_role(r);
            role_to_phys.erase(r);
        }
    }
#ifdef PROTEAN_PERF
    fp_t t = timer.clk_end();
    std::cout << "[ rcr ] role merging took " << t*1e-9 << "s" << std::endl;
#endif
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

#ifdef PROTEAN_PERF
    timer.clk_start();
#endif

    uptr<Graph<int_v_t, int_e_t>> interaction_graph = std::make_unique<Graph<int_v_t, int_e_t>>();
    size_t _id = 0;
    for (sptr<tanner::vertex_t> tpq : tanner_graph->get_checks()) {
        sptr<raw_vertex_t> rpq = raw_connection_network->v_tanner_raw_map.at(tpq);
        // Make interaction graph vertex.
        sptr<int_v_t> iv = interaction_graph->make_and_add_vertex(_id++);
        iv->check = rpq;
        // Get all qubits involved in check.
        RawNetwork::parity_support_t& supp = raw_connection_network->get_support(rpq);
        iv->support = supp.all;
    }
    // Now, compute interaction graph edges.
    for (size_t i = 0; i < _id; i++) {
        sptr<int_v_t> iv = interaction_graph->get_vertex(i);
        for (size_t j = i+1; j < _id; j++) {
            sptr<int_v_t> iw = interaction_graph->get_vertex(j);

            sptr<int_e_t> ie = interaction_graph->make_edge(iv, iw);
            for (sptr<raw_vertex_t> rx : iv->support) {
                if (rx->qubit_type == raw_vertex_t::type::data) continue;

                sptr<phys_vertex_t> px = role_to_phys.at(rx);
                for (sptr<raw_vertex_t> ry : iw->support) {
                    if (rx == ry && !rx->is_check()) continue;
                    sptr<phys_vertex_t> py = role_to_phys.at(ry);
                    if (px == py) ie->conflicts.emplace_back(rx, ry, px);
                }
            }
            if (ie->conflicts.size()) interaction_graph->add_edge(ie);
        }
    }
    // Clear all the cycle role maps for each physical qubit.
    for (sptr<phys_vertex_t> pv : get_vertices()) pv->clear_roles();
    // Color the interaction graph.
    std::map<sptr<int_v_t>, int> color_map;
    size_t cycles = static_cast<size_t>(k_coloring_rlf(interaction_graph.get(), color_map));
    for (size_t c = 0; c <= cycles; c++) {
        for (sptr<int_v_t> iv : interaction_graph->get_vertices()) {
            if (color_map[iv] != c) continue;
            for (sptr<raw_vertex_t> rv : iv->support) {
                if (!role_to_phys.count(rv) || role_to_phys.at(rv) == nullptr) {
                    std::cerr << "Role does not exist: " << print_v(rv) << "\n"
                        << "\tcontaining check: " << print_v(iv->check) << "\n"
                        << "\traw_net contains role: " << raw_connection_network->contains(rv) << std::endl;
                    exit(1);
                }
                sptr<phys_vertex_t> pv = role_to_phys.at(rv);
                if (!pv->role_set.count(rv)) {
                    if (rv->qubit_type == raw_vertex_t::type::data) {
                        pv->push_back_role(rv, 0);
                    } else {
                        pv->push_back_role(rv, c);
                    }
                }
            }
        }
    }
#ifdef PROTEAN_PERF
    t = timer.clk_end();
    std::cout << "[ rcr ] interference analysis took " << t*1e-9 << "s" << std::endl;
#endif
    return true;
}

void
PhysicalNetwork::finalize() {
    // Run the scheduler.
    Scheduler sch(this);
    x_memory = sch.run(config.rounds, true);
    z_memory = sch.run(config.rounds, false);
    round_latency = sch.get_depth_as_time();
    round_cnots = sch.get_depth_as_cx_opcount();
    // Run the yield simulator.
    /*
    YieldSimulator ysim(this);
    std::vector<fp_t> freq_list;
    for (fp_t f = config.min_qubit_frequency; 
            f <= config.max_qubit_frequency+0.1; 
            f += config.qubit_frequency_step) 
    {
        freq_list.push_back(f);
    }
    ysim.assign(config.fabrication_precision, freq_list);
    ysim.assign(config.fabrication_precision, freq_list);
    expected_collisions = ysim.est_mean_collisions(config.fabrication_precision, 10000);
    */
}

}   // protean
}   // qontra
