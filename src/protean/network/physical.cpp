/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "protean/network.h"

#include <lemon/matching.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

PhysicalNetwork(TannerGraph& tgr, uint max_conn)
    :Graph(),
    raw_connection_network(tgr),
    role_to_phys(),
    // Planarity tracking:
    planar_repr(),
    v_phys_to_lemon(),
    v_lemon_to_phys(),
    e_phys_to_lemon(),
    e_lemon_to_phys(),
    // Other tracking:
    occupied_tsvs(),
    bulk_degree_map(),
    // Other variables:
    max_connectivity(max_conn)
{
    // Create a corresponding physical qubit for every vertex in raw_connection_network.
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        sptr<phys_vertex_t> pv = std::make_shared<>();
        pv->id = rv->id;
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
PhysicalNetwork::make_flags(bool for_x_parity) {
    using namespace lemon;
    // Set constants based on for_x_parity.
    const raw_vertex_t::type desired_role_type =
        for_x_parity ? raw_vertex_t::type::xparity : raw_vertex_t::type::zparity;

    ListGraph matching_graph;
    ListGraph::EdgeMap<int> edge_weights;

    BijectiveMap<sptr<phys_vertex_t>, ListGraph::Node> data_lemon_map;

    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        // We only care about data qubits.
        if (rv->qubit_type != raw_vertex_t::type::data) continue;
        sptr<phys_vertex_t> pv = role_to_phys[rv];
        if (!data_lemon_map.count(pv)) {
            data_lemon_map.put(pv, matching_graph.addNode());
        }
    }
    // Now for each vertex, make an edge with other vertices that share common neighbors.
    TwoLevelMap<sptr<phys_vertex_t>, sptr<phys_vertex_t>, std::set<sptr<phys_vertex_t>>>
        common_neighbors_map;

    std::set<sptr<phys_vertex_t>> visited;
    for (auto& p1 : data_lemon_map) {
        sptr<phys_vertex_t> pv = p1.first;
        auto vnode = p1.second;

        // We care about the number of common neighbors of a given parity type (x or z, as
        // specified in the function input).
        std::set<sptr<phys_vertex_t>> p_neighbors;
        for (sptr<phys_vertex_t> x : get_neighbors(pv)) {
            if (x->has_role_of_type(desired_role_type)) p_neighbors.insert(x);
        }

        for (auto& p2 : data_lemon_map) {
            if (p1 == p2) continue;
            sptr<phys_vertex_t> pw = p2.first;
            auto wnode = p2.second;
            // Do not double count (check if pw is visited already).
            if (visited.count(pw)) continue;
            // Get common neighbors:
            std::set<sptr<phys_vertex_t>> common;
            for (sptr<phys_vertex_t> x : get_neighbors(pw)) {
                if (p_neighbors.count(x)) common.insert(x);
            }
            // Create edge if common is not empty.
            if (!common.empty()) {
                auto edge = ListGraph::addEdge(vnode, wnode);
                edge_weights[edge] = static_cast<int>(common.size());

                tlm_put(common_neighbors_map, pv, pw, common);
                tlm_put(common_neighbors_map, pw, pv, common);
            }
        }
        visited.count(pv);
    }
    // Now, we can perform Max-Weight Matching.
    MaxWeightPerfectMatching maxwpm(matching_graph, edge_weights);
    maxwpm.run();
    for (auto& pair : data_lemon_map) {
        sptr<phys_vertex_t> pv = p1.first;
        auto vnode = p1.second;
        auto wnode = maxwpm.mate(vnode);
        sptr<phys_vertex_t> pw = data_lemon_map[wnode];
        // So flag creation will happen as follows:
        //
        //  (1) For every common neighbor between pv and pw, we will get all parity qubit
        //  roles and create a flag. Note that even if we did X parity qubits right now,
        //  if an adjacent parity qubit has a role as a Z parity qubit, we will attach the
        //  flag to that Z parity qubit as well (this is space efficient). Note that these
        //  flags are in the raw_connection_network. This also happens for each role in 
        //  px and py (provided that the roles rx and ry share an edge with the parity role).
        //
        //  (2) Each raw flag will have a single physical flag, which is connected to pv
        //  and pw and all qubits in common.
        sptr<phys_vertex_t> pfq = std::make_shared<>();
        std::vector<sptr<phys_vertex_t>> pfq_make_edge_with_list{pv, pw};

        auto common = common_neighbors_map[pv][pw];
        for (sptr<phys_vertex_t> ppx : common) {
            for (sptr<raw_vertex_t> rv : pv->role_set) {
                for (sptr<raw_vertex_t> rw : pw->role_set) {
                    for (sptr<raw_vertex_t> rpx : ppx->role_set) {
                        // Not possible to make a flag if the following is not true:
                        if (!raw_connection_network.has_edge(rv, rpx)
                            || !raw_connection_network.has_edge(rw, rpx))
                        {
                            continue;
                        }
                        sptr<raw_vertex_t> rfq = raw_connection_network.add_flag(rv, rw, rpx);
                        // The cycle of rfq is the max of rv, rw, and rpx cycles.
                        size_t rv_cycle = pv->cycle_role_map.at(rv),
                                rw_cycle = pw->cycle_role_map.at(rw),
                                rpx_cycle = ppx->cycle_role_map.at(rpx);
                        size_t rfq_cycle = std::max(rv_cycle, rw_cycle, rpx_cycle);
                        pfq->add_role(rfq, rfq_cycle);
                    }
                }
            }
            // Delete edge between pv (pw) and ppx.
            sptr<phys_edge_t> e1 = get_edge(pv, ppx);
            sptr<phys_edge_t> e2 = get_edge(pw, ppx);
            delete_edge(e1);
            delete_edge(e2);

            pfq_make_edge_with_list.push_back(ppx);
        }
        // Now, connect pfq to pv, pw, and all in common.
        for (sptr<phys_vertex_t> x : pfq_make_edge_with_list) {
            sptr<phys_edge_t> e = std::make_shared<>();
            e->src = pfq;
            e->dst = x;
            e->is_undirected = true;
            add_edge(e);
        }
    }
}

}   // protean
}   // qontra
