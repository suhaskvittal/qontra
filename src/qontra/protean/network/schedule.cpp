/*
 *  author: Suhas Vittal
 *  date:   7 January 2024
 * */

#include "qontra/protean/network.h"

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

template <class V> using v_pair_t=std::pair<sptr<V>, sptr<V>>;

inline v_pair_t<raw_vertex_t>
orient(sptr<raw_vertex_t> v, sptr<raw_vertex_t> w) {
    if (v->qubit_type == raw_vertex_t::type::xparity) {
        return std::make_pair(v, w);
    } else {
        return std::make_pair(w, v);
    }
}

inline v_pair_t<raw_vertex_t>
orient_relative(sptr<raw_vertex_t> v, sptr<raw_vertex_t> w, sptr<raw_vertex_t> to) {
    if (to->qubit_type == raw_vertex_t::type::xparity) {
        return std::make_pair(v, w);
    } else {
        return std::make_pair(w, v);
    }
}


qes::Program<>
PhysicalNetwork::make_schedule(size_t rounds) {
    qes::Program<> program;
    // Recompute cycle role maps.
    recompute_cycle_role_maps();
    // The goal of this algorithm is to perform a CNOT
    // for every edge in raw_connection_network. The steps are:
    //  (1) Assign roles to each parity role.
    //  (2) For each cycle, we identify all edges that can occur in that cycle
    //      (cycles of both endpoints must be at least the current cycle). Then:
    //      (a) Preparation: do all flag edges.
    //          --> Propagate CNOTs from Z flags and X checks.
    //      (b) CNOTs: do all data edges
    //          --> Propagate CNOTs from data qubits to Z flags/checks and X flags/checks to CNOTs.
    //      (c) Teardown: do all flag edges again.
    //          --> Propagate CNOTs as before.
    //      (d) Measure + Reset + Detection Events. Cycle ends now.
    //      
    //      We must track propagation as if we cannot do an edge in the current cycle, we must push
    //      back operations to the next cycle.
    struct parity_support_t {
        std::set<sptr<raw_vertex_t>> data;
        std::set<sptr<raw_vertex_t>> flags;
    };
    std::set<sptr<raw_vertex_t>> all_checks;
    std::map<sptr<raw_vertex_t>, parity_support_t> parity_support_map;
    // Populate the parity support map.
    TannerGraph& tanner_graph = raw_connection_network.tanner_graph;
    for (sptr<tanner::vertex_t> tpq : tanner_graph.get_checks()) {
        sptr<raw_vertex_t> rpq = raw_connection_network.v_tanner_raw_map.at(tpq);
        all_checks.insert(rpq);

        parity_support_t support;
        for (sptr<tanner::vertex_t> tdq : tanner_graph.get_neighbors(tpq)) {
            sptr<raw_vertex_t> rdq = raw_connection_network.v_tanner_raw_map.at(tdq);
            support.data.insert(rdq);
        }
        support.flags.insert(raw_connection_network.flag_ownership_map[rpq].begin(),
                            raw_connection_network.flag_ownership_map[rpq].end());
        parity_support_map[rpq] = support;
    }
    // Now, let's make that schedule.
    //
    // Data structures:
    //  (1) stage_map -- track where each check is in syndrome extraction
    //                      0 = preparation, 1 = CX, 2 = teardown, 3 = m+r+e
    //  (2) visited_edge_map -- track which role edges have been used and in which stage. If
    //                      an edge has already been used in a stage, we know not to repeat it.
    //  (3) h_gate_stage_map -- track if we have done an H gate for a stage.
    std::map<sptr<raw_vertex_t>, int>       stage_map;
    std::map<sptr<raw_edge_t>, int>         visited_edge_map;
    std::map<sptr<raw_vertex_t>, int>       h_gate_stage_map;
    
    // Set visited_edge_map and h_gate_stage_map to all -1
    for (auto rv : raw_connection_network.get_vertices()) h_gate_stage_map[rv] = -1;
    for (auto re : raw_connection_network.get_edges()) visited_edge_map[re] = -1;

    vtils::TwoLevelMap<sptr<raw_vertex_t>, sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>>
        proxy_walk_memo_map;

    size_t curr_cycle = 0;

    bool done;
    while (!done) {
        Program<> preparation, body;
        // Preparation:
        // 
        // We want to perform an H gate on all Z flags and X checks, along with a CNOT
        // from the Z flags to Z checks and X checks to X flags.
        {
            std::vector<uint64_t> h_operands;
            for (auto& p : parity_support_map) {
                sptr<raw_vertex_t> rpq = p.first;
                parity_support_t& support = p.second;
                if (stage_map[rpq] != 0) continue;
                // Attempt to perform any H gates.
                std::vector<sptr<raw_vertex_t>> _h_operands;
                if (rpq->qubit_type == raw_vertex_t::type::xparity) {
                    _h_operands.push_back(rpq);
                } else {
                    _h_operands.insert(_h_operands.end(),
                            parity_support_map[rpq].flags.begin(),
                            parity_support_map[rpq].flags.end());
                }
                // Perform the H gates if possible.
                for (sptr<raw_vertex_t> rx : _h_operands) {
                    if (h_gate_stage_map[rx] == 0) continue; // already done.
                    sptr<phys_vertex_t> px = role_to_phys[rx];
                    // Check if it is time for this qubit.
                    if (px->cycle_role_map[rx] > curr_cycle) continue;
                    h_operands.push_back(px->id);
                    h_gate_stage_map[rx] = 0;
                }
            }
            preparation.emplace_back("h", h_operands);
            preparation.back().put("timing_error");
            // Now, we must perform the CNOTS.
            std::map<sptr<raw_vertex_t>, size_t> k_map;
            while (1) {
                std::vector<uint64_t> cx_operands;
                for (auto& p : parity_support_map) {
                    sptr<raw_vertex_t> rpq = p.first;
                    parity_support_t& support = p.second;
                    if (stage_map[rpq] != 0) continue;
                    if (k_map[rpq] >= support.flags.size()) continue;
                    // Attempt to perform any CX gates.
                    v_pair_t<raw_vertex_t> cx_edge = std::make_pair(nullptr, nullptr);
                    sptr<raw_vertex_t> rfq = *(support.flags.begin() + k_map[rpq]);
                    // Check if rpq and rfq share an edge. If so, all good.
                    if (raw_connection_network.has_edge(rpq, rfq)) {
                        sptr<raw_edge_t> re = raw_connection_network.get_edge(rpq, rfq);
                        if (visited_edge_map[re] != 0) {
                            cx_edge = orient(rpq, rfq);
                            k_map[rpq]++;
                        }
                        
                    } else {
                        // Now, here's where it becomes tricky. rpq and rfq are separated
                        // by proxies. Now, we have to search for the proxy.
                        if (!proxy_memo_map.count(rpq) && !proxy_memo_map[rpq].count(rfq)) {
                            // We must search.
                            for (sptr<raw_vertex_t> rprx : raw_connection_network.get_neighbors(rpq)) {
                                if (rprx->qubit_type != raw_vertex_t::type::proxy) continue;
                                std::vector<sptr<raw_vertex_t>> walk_path;
                                if (rfq == raw_connection_network.proxy_walk(rpq, rprx, walk_path)) {
                                    // We found the proxy qubit.
                                    vtils::tlm_put(proxy_memo_map, rpq, rfq, walk_path);
                                    vtils::tlm_put(proxy_memo_map, rfq, rpq, walk_path);
                                }
                            }
                            // We will should not have to do this again.
                        }
                        std::vector<sptr<raw_vertex_t>> path = proxy_memo_map[rpq][rfq];
                        size_t i;
                        for ( ; i < path.size(); i++) {
                            sptr<raw_edge_t> re = raw_connection_network.get_edge(path[i-1], path[i]);
                            if (visited_edge_map[re] != 0) {
                                cx_edge = orient_relative(path[i-1], path[i], rpq);
                                break;
                            }
                        }
                        if (i == path.size()-1) k_map[rpq]++
                    }
                    if (cx_edge.first == nullptr) {
                        // Only way this happens is if there is no edge to do, so increment k_map.
                        k_map[rpq]++;
                        continue;
                    };
                    // Otherwise, proceed and check if we can perform the edge.
                    sptr<raw_vertex_t> rx = cx_edge.first,
                                        ry = cx_edge.second;
                    sptr<phys_vertex_t> px = role_to_phys[rx],
                                        py = role_to_phys[ry];
                    if (px->cycle_role_map.at(rx) > curr_cycle
                        || py->cycle_role_map.at(ry) > curr_cycle)
                    {
                        continue;
                    }
                    sptr<raw_edge_t> re = raw_connection_network.get_edge(rx, ry);
                    visited_edge_map[re] = 0;
                    // Perform the cx.
                    cx_operands.push_back(px->id, py->id);
                }
                if (cx_operands.empty()) break;
                preparation.emplace_back("cx", cx_operands);
            }
            // Now, we are done with preparation. We set stage_map to 1 if:
            //  (1) h_gate_stage_map of Z flags/X checks are set to 0, and
            //  (2) k_map is set to # of flags.
            for (auto& p : parity_support_map) {
                sptr<raw_vertex_t> rpq = p.first;
                parity_support_t& support = p.second;
                if (stage_map[rpq] != 0) continue;
                // Check condition 1.
                bool cond_1 = true;
                for (sptr<raw_vertex_t> rfq : support.flags) {
                    cond_1 &= (h_gate_stage_map[rfq] == 0);
                }
                // Condition 2 is a one liner.
                bool cond_2 = (k_map[rpq] == support.flags.size());
                if (cond_1 && cond_2) {
                    stage_map[rpq] = 1;
                }
            }
        }
        // CX Body
        {

        }
    }
}

}   // protean
}   // qontra
