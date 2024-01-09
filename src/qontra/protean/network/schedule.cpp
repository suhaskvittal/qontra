/*
 *  author: Suhas Vittal
 *  date:   7 January 2024
 * */

#include "qontra/protean/network.h"

#include <vtils/linprog/manager.h>
#include <vtils/two_level_map.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

template <class V> using v_pair_t=std::pair<sptr<V>, sptr<V>>;

static const v_pair_t<raw_vertex_t> NULL_CX_EDGE = std::make_pair(nullptr, nullptr);

Scheduler::Scheduler(PhysicalNetwork* n_p)
    :all_checks(),
    parity_support_map(),
    stage_map(),
    visited_edge_map(),
    h_gate_stage_map(),
    proxy_walk_memo_map(),
    net_p(n_p),
    cycle(0)
{
    // Initialize all structures.
    TannerGraph& tanner_graph = net_p->raw_connection_network.tanner_graph;
    // First, all_checks and parity_support_map
    for (sptr<tanner::vertex_t> tpq : tanner_graph.get_checks()) {
        sptr<raw_vertex_t> rpq = raw_connection_network.v_tanner_raw_map.at(tpq);
        all_checks.insert(rpq);

        parity_support_t support;
        for (sptr<tanner::vertex_t> tdq : tanner_graph.get_neighbors(tpq)) {
            sptr<raw_vertex_t> rdq = raw_connection_network.v_tanner_raw_map.at(tdq);
            support.data.insert(rdq);
        }
        support.flags.insert(net_p->raw_connection_network.flag_ownership_map[rpq].begin(),
                            net_p->raw_connection_network.flag_ownership_map[rpq].end());
        parity_support_map[rpq] = support;
    }
    // Now, the tracking structures. Just set default values.
    //
    // Set visited_edge_map and h_gate_stage_map to all -1
    for (auto rv : net_p->raw_connection_network.get_vertices()) h_gate_stage_map[rv] = -1;
    for (auto re : net_p->raw_connection_network.get_edges()) visited_edge_map[re] = -1;
    // No need to initialize stage_map explicitly, default value is 0.
}

qes::Program<>
Scheduler::build_preparation() {
    qes::Program<> program;

    std::set<sptr<raw_vertex_t>> checks_this_stage = get_checks_at_stage(PREP_STAGE);
    // First, perform the H gates.
    std::vector<uint64_t> h_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        h_operands.insert_range(preparation_get_h_operands(rpq));
    }
    program.emplace_back("h", h_operands);
    program.back().put("timing_error");
    // Now, we must perform the CNOTS.
    std::map<sptr<raw_vertex_t>, size_t> k_map;
    while (1) {
        std::vector<uint64_t> cx_operands;
        for (sptr<raw_vertex_t> rpq : checks_this_stage) {
            cx_operands.insert_range(preparation_get_cx_operands(rpq, k_map[rpq]));
        }
        if (cx_operands.empty()) break;
        program.emplace_back("cx", cx_operands);
    }
    // Now, we are done with preparation. We set stage_map to 1 if:
    //  (1) h_gate_stage_map of Z flags/X checks are set to 0, and
    //  (2) k_map is set to # of flags.
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        sptr<raw_vertex_t> rpq = p.first;
        auto support = parity_support_map[rpq];
        // Check condition 1.
        bool cond_1 = true;
        for (sptr<raw_vertex_t> rfq : support.flags) {
            cond_1 &= (h_gate_stage_map[rfq] == 0);
        }
        // Condition 2 is a one liner.
        bool cond_2 = (k_map[rpq] == support.flags.size());
        if (cond_1 && cond_2) {
            stage_map[rpq] = BODY_STAGE;
        }
    }
    return program;
}

qes::Program<>
Scheduler::build_body() {
    using namespace vtils;

    std::set<sptr<raw_vertex_t>> checks_this_stage = get_checks_at_stage(BODY_STAGE);
    // Here, we will build a linear program (a MIP specifically) to compute the schedule
    // based off the partial supports.
    
    // Maps parity qubit -> check schedule array (nullptr is where a data qubit CNOT is not scheduled).
    std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>>
        data_cnot_schedules;
    // Maps data qubit -> vector of (time, originating_check).
    std::map<sptr<raw_vertex_t>, std::vector<std::pair<size_t, sptr<raw_vertex_t>>>> 
        existing_data_cnot_times;

    // First, we need the max check operator weight for this stage.
    size_t check_operator_max_weight = 0;
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        check_operator_max_weight = std::max(parity_support_map[rpq].data.size(), check_operator_max_weight);
    }
    const size_t upper_bound = 2*check_operator_max_weight;
    // Now, we can make the schedules.
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        bool is_x_check = (rpq->qubit_type == raw_vertex_t::type::xparity);
        auto partial_support = body_get_partial_data_support(rpq);
        // Build the LP.
        CPXLPManager<sptr<raw_vertex_t>> LP;
        lp_var_t max_of_all = LP.add_slack_var(0, upper_bound, lp_var_t::bounds::both, lp_var_t::domain::integer);
        for (sptr<raw_vertex_t> rdq : partial_support) {
            // Build variable.
            lp_var_t x = LP.add_var(rdq, 0, upper_bound, lp_var_t::bounds::both, lp_var_t::domain::integer);

            lp_constr_t max_con(max_of_all - x, 0.0, lp_constr_t::direction::ge);
            LP.add_constraint(max_con);
        }
        // Add uniqueness constraints
        for (sptr<raw_vertex_t> rq1 : partial_support) {
            lp_var_t x = LP.get_var(rq1);
            for (sptr<raw_vertex_t> rq2 : partial_support) {
                if (rq1 == rq2) continue;
                lp_var_t y = LP.get_var(rq2);
                lp_constr_t con(x, y, lp_constr_t::direction::neq);
                LP.add_constraint(con);
            }
        }
        // Add stabilizer constraints.
        std::map<sptr<raw_vertex_t>, lp_expr_t> ind_sum_map;
        for (sptr<raw_vertex_t> rdq : partial_support) {
            lp_var_t x = LP.get_var(rdq);
            for (auto& p : existing_data_cnot_times[rdq]) {
                size_t t = p.first;
                sptr<raw_vertex_t> rx = p.second;
                const bool is_from_x = (rx->qubit_type == raw_vertex_t::type::xparity);
                // Add time-conflict constraint.
                lp_constr_t con1(x, t, lp_constr_t::direction::neq);
                LP.add_constraint(con1);
                // if is_from_x != is_x_check, then add commutativity constraints.
                if (is_from_x == is_x_check) continue;
                lp_expr_t& ind_sum = ind_sum_map[rx];
                
                const double M = 100000;
                lp_var_t y = LP.add_slack_var(0, 1, lp_var_t::bounds::both, lp_var_t::domain::binary);
                // x - t <= My and t - x <= M(1-y)
                lp_constr_t con2(x - t, M*y, lp_constr_t::direction::le);
                lp_constr_t con3(x - t, M - M*y, lp_constr_t::direction::le);
                LP.add_constraint(con2);
                LP.add_constraint(con3);
                ind_sum += y;
            }
        }
        for (auto p : ind_sum_map) {
            lp_expr_t ind_sum = p.second;
            lp_var_t y = LP.add_slack_var(0, 0, lp_var_t::bounds::lower, lp_var_t::domain::integer);
            lp_constr_t con(ind_sum, 2*y, lp_constr_t::direction::eq);
            LP.add_constraint(con);
        }
        // Solve the LP.
        LP.build(max_of_all, false);
        double obj;
        int solstat;
        int status;
        if ((status=LP.solve(&obj, &solstat))) {
            std::cerr << "build_body: program is infeasible, solstat = " << solstat
                << ", status = " << status << std::endl;
            exit(1);
        }
        // Update the schedule.
        std::vector<sptr<raw_vertex_t>> check_sch(static_cast<size_t>(obj), nullptr);
        for (sptr<raw_vertex_t> rdq : partial_support) {
            size_t t = LP.get_value(rdq);
            existing_data_cnot_times[rdq].push_back(std::make_pair(t, rpq));
            check_sch[t] = rdq;
        }
        data_cnot_schedules[rpq] = check_sch;
    }
    // Now, that we have the CNOT schedules of all the partial supports, we need to execute all possible
    // CNOTs for this cycle.

}

std::vector<sptr<raw_vertex_t>>
Scheduler::get_proxy_walk_path(sptr<raw_vertex_t> src, sptr<raw_vertex_t> dst) {
    static vtils::TwoLevelMap<sptr<raw_vertex_t>, sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>>
        proxy_walk_memo_map;

    if (!proxy_memo_map.count(src) && !proxy_memo_map[dst].count(rfq)) {
        // We must search.
        for (sptr<raw_vertex_t> rprx : net_p->raw_connection_network.get_neighbors(src)) {
            if (rprx->qubit_type != raw_vertex_t::type::proxy) continue;
            std::vector<sptr<raw_vertex_t>> walk_path;
            if (dst == net_p->raw_connection_network.proxy_walk(src, dst, walk_path)) {
                // We found the proxy qubit.
                vtils::tlm_put(proxy_memo_map, src, dst, walk_path);
                std::reverse(walk_path.begin(), walk_path.end());
                vtils::tlm_put(proxy_memo_map, dst, src, walk_path);
            }
        }
    }
    return proxy_walk_memo_map[src][dst];
}

sptr<raw_vertex_t>
Scheduler::test_and_get_other_endpoint_if_ready(
        sptr<raw_vertex_t> endpoint,
        std::vector<sptr<raw_vertex_t>> path) 
{
    for (size_t i = 1; i < path.size(); i++) {
        sptr<raw_vertex_t> rx = path[i-1],
                            ry = path[i];
        sptr<raw_edge_t> re = raw_net.get_edge(rx, ry);
        if (rx == endpoint || ry == endpoint) {
            if (!is_good_for_current_cycle(rx) || !is_good_for_current_cycle(ry)) return nullptr;
            if (visited_edge_map[re] < BODY_STAGE) {
                // Then the edge is not done.
                return rx == endpoint ? ry : rx;
            }
        } else {
            if (visited_edge_map[re] < BODY_STAGE) {
                return nullptr;
            }
        }
    }
    return nullptr;
}

std::vector<uint64_t>
Scheduler::preparation_get_h_operands(sptr<net::raw_vertex_t> rpq) {
    std::vector<uint64_t> operands;

    parity_support_t& support = parity_support_map[rpq];
    // Collect potential operands.
    if (rpq->qubit_type == raw_vertex_t::type::xparity) {
        operands.push_back(rpq);
    } else {
        operands.insert(operands.end(), support.flags.begin(), support.flags.end());
    }
    // Perform the H gates if possible.
    for (auto it = operands.begin(); it != operands.end(); ) {
        sptr<raw_vertex_t> rx = *it;
        if (h_gate_stage_map[rx] == PREP_STAGE) {
            // We have already done this gate in this stage.
            it = operands.erase(it);
        } else {
            sptr<phys_vertex_t> px = net_p->role_to_phys[rx];
            if (!is_good_for_curr_cycle(rx)) {
                // We cannot do this operation yet.
                it = operands.erase(it);
            } else {
                operands.push_back(px->id);
                h_gate_stage_map[rx] = PREP_STAGE;
                it++;
            }
        }
    }
    return operands;
}

std::vector<uint64_t>
PhysicalNetwork::preparation_get_cx_operands(sptr<net::raw_vertex_t> rpq, size_t& k) {
    // k is corresponds to the operand index.
    std::vector<uint64_t> operands;

    parity_support_t& support = parity_support_map;
    // Check if k is valid.
    if (k >= support.flags.size()) return operands;
    // Now, attempt to get CX operands.
    v_pair_t<raw_vertex_t> cx_edge = NULL_CX_EDGE;
    // Get k-th flag qubit.
    sptr<raw_vertex_t> rfq = *(support.flags.begin() + k);
    // There are two possiblities.
    //  (1) rpq and rfq are directly connected. Then, the CX is straightforward.
    //  (2) rpq and rfq are not directly connected. We need to do a proxy walk.
    if (net_p->raw_connection_network.has_edge(rpq, rfq)) {
        sptr<raw_edge_t> re = net_p->raw_connection_network.get_edge(rpq, rfq);
        if (visited_edge_map[re] < PREP_STAGE) {
            cx_edge = orient(rpq, rfq);
            k++;
        }
    } else {
        // If rpq is a Z check, then the CNOTs run from rfq to rpq.
        std::vector<sptr<raw_vertex_t>> proxy_path = rpq->qubit_type == raw_vertex_t::type::xparity 
                                                        ? get_proxy_walk_path(rpq, rfq)
                                                        : get_proxy_walk_path(rfq, rpq);
        for (size_t i = 0; i < path.size(); i++) {
            sptr<raw_edge_t> re = net_p->raw_connection_network.get_edge(path[i-1], path[i]);
            if (visited_edge_map[re] < PREP_STAGE) {
                // No need to orient.
                cx_edge = std::make_pair(path[i-1], path[i]);
                // If we are at the end of the path, then we are done with this flag, so move on (inc k).
                if (i == path.size() - 1) k++;
                break;
            }
        }
    }
    // This occurs when the edge is already done. Here, we such call this function recursively.
    if (cx_edge == NULL_CX_EDGE) return preparation_get_cx_operands(rpq, ++k);
    sptr<raw_vertex_t> rx = cx_edge.first,
                       ry = cx_edge.second;
    // Check if the edge is doable.
    if (!is_good_for_curr_cycle(rx) || !is_good_for_curr_cycle(ry)) {
        return operands;
    }
    // Update the visited_edge_map.
    sptr<raw_edge_t> re = net_p->raw_connection_network.get_edge(rx, ry);
    visited_edge_map[re] = PREP_STAGE;
    // Return the operands.
    sptr<phys_vertex_t> px = net_p->role_to_phys[rx],
                        py = net_p->role_to_phys[ry];
    operands.push_back(px->id);
    operands.push_back(py->id);
    return operands;
}

std::set<sptr<raw_vertex_t>>
Scheduler::body_get_partial_data_support(sptr<raw_vertex_t> rpq) {
    std::set<sptr<raw_vertex_t>> partial_support;

    RawNetwork& raw_net = net_p->raw_connection_network;

    parity_support_t& support = parity_support_map[rpq];
    for (sptr<raw_vertex_t> rdq : support.data) {
        // Four possibilites:
        //  (1) rdq is directly connected to rpq (no flag)
        //  (2) rdq is indirectly connected to rpq (proxy walk)
        //  (3) rdq is directly connected to a flag
        //  (4) rdq is indirectly connected to a flag
        sptr<raw_vertex_t> other = nullptr;
        if (raw_net.has_edge(rdq, rpq)) {
            // Situation 1
            other = rpq;
        } else if (raw_net.flag_assignment_map[rpq].count(rdq)) {
            sptr<raw_vertex_t> rfq = raw_net[rpq][rdq];
            // Check if rdq is connected to rfq.
            if (raw_net.has_edge(rdq, rfq)) {
                // Situation 3
                other = rfq;
            } else {
                // Situation 4.
                auto proxy_walk_path = rpq->qubit_type == raw_vertex_t::type::xparity
                                        ? get_proxy_walk_path(rfq, rdq)
                                        : get_proxy_walk_path(rdq, rfq);
                other = test_and_get_edge_if_ready(rdq, proxy_walk_path);
                // Now, we need to see if rdq's edge can be done (or if it is already done).
            }
        } else {
            // Situation 2.
            auto proxy_walk_path = rpq->qubit_type == raw_vertex_t::type::xparity
                                    ? get_proxy_walk_path(rpq, rdq)
                                    : get_proxy_walk_path(rdq, rpq);
            other = test_and_get_edge_if_ready(rdq, proxy_walk_path);
        }
        // If other == nullptr, this means that we are not ready to do the edge in the proxy walk.
        if (other == nulltpr) continue;
        if (!is_good_for_curr_cycle(rdq) || !is_good_for_curr_cycle(other)) continue;
        sptr<raw_edge_t> re = raw_net.get_edge(rdq, other);
        if (visited_edge_map[re] < BODY_STAGE) {
            partial_support.insert(rdq);
        }
    }
    return partial_support;
}

std::vector<uint64_t>
Scheduler::body_get_cx_operands(
        sptr<raw_vertex_t> rpq,
        const std::vector<sptr<raw_vertex_t>>& cnot_sch,
        size_t& k)
{
    std::vector<uint64_t> operands;
    if (k >= cnot_sch.size()) return operands;

    sptr<raw_vertex_t> rdq = cnot_sch.at(k);
    // Check if there is a data qubit at time step k.
    if (rdq != nullptr) {
        // Then we are obligated to perform this CNOT.
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
}

}   // protean
}   // qontra
