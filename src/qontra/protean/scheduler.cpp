/*
 *  author: Suhas Vittal
 *  date:   7 January 2024
 * */

#include "qontra/protean/scheduler.h"

#include <vtils/linprog/manager.h>
#include <vtils/two_level_map.h>
#include <vtils/utility.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;
using namespace vtils;

template <class V> using v_pair_t=std::pair<sptr<V>, sptr<V>>;
static const v_pair_t<raw_vertex_t> NULL_CX_EDGE = std::make_pair(nullptr, nullptr);

inline void
safe_emplace_back(qes::Program<>& program, std::string name, std::vector<uint64_t> operands) {
    if (operands.empty()) return;
    program.emplace_back(name, operands);
}

Scheduler::Scheduler(PhysicalNetwork* n_p)
    :all_checks(),
    stage_map(),
    visited_edge_map(),
    h_gate_stage_map(),
    flag_stage_map(),
    proxy_reset_map(),
    proxy_occupied_map(),
    active_role_map(),
    remaining_roles_map(),
    data_stage_map(),
    data_cx_map(),
    cx_in_use_set(),
    meas_ctr_map(),
    meas_ctr(0),
    running_ind_sum_map(),
    running_ind_ctr_map(),
    cycle(0),
    cx_return_status(0),
    net_p(n_p)
{
    // Initialize all structures.
    RawNetwork& raw_net = net_p->raw_connection_network;
    TannerGraph& tanner_graph = raw_net.tanner_graph;
    // First, all_checks
    for (sptr<tanner::vertex_t> tpq : tanner_graph.get_checks()) {
        sptr<raw_vertex_t> rpq = raw_net.v_tanner_raw_map.at(tpq);
        all_checks.insert(rpq);
    }
    // Now, the tracking structures. Just set default values.
    for (auto rv : raw_net.get_vertices()) {
        h_gate_stage_map[rv] = -1;
        flag_stage_map[rv] = -1;
    }
    for (auto rpq : all_checks) {
        for (auto rdq : get_support(rpq).data) {
            tlm_put(data_stage_map, rpq, rdq, -1);
        }
    }
    for (auto re : raw_net.get_edges()) visited_edge_map[re] = -1;

    for (auto pv : net_p->get_vertices()) {
        std::vector<sptr<raw_vertex_t>> roles(pv->role_set.begin(), pv->role_set.end());
        std::sort(roles.begin(), roles.end(),
                [&] (auto rx, auto ry)
                {
                    // We want the lower cycle roles towards the back.
                    return pv->cycle_role_map.at(rx) > pv->cycle_role_map.at(ry);
                });
        remaining_roles_map[pv] = std::move(roles);
    }
    // No need to initialize stage_map explicitly, default value is 0.
}

qes::Program<>
Scheduler::run() {
    qes::Program<> program;
    
    bool done;
    do {
        size_t first_inst_of_cycle = program.size();
        build_preparation(program);
        build_body(program);
        build_teardown(program);
        build_measurement(program);
        
        cycle++;
        // Check if we are done (all checks are in DONE_STAGE).
        done = true;
        for (sptr<raw_vertex_t> rpq : all_checks) {
            done &= (stage_map[rpq] == DONE_STAGE);
        }
        // Add timing_error annotation
        if (first_inst_of_cycle < program.size()) {
            program[first_inst_of_cycle].put("timing_error");
        }
#ifdef PROTEAN_DEBUG
        // Print status of each parity check:
        std::cerr << "[ scheduler ] parity check status:";
        for (sptr<raw_vertex_t> rpq : all_checks) {
            std::cerr << "\n" << print_v(rpq) << " --> " << stage_map.at(rpq);
        }
        std::cerr << std::endl;
#endif
    } while (!done);
    return program;
}

void
Scheduler::build_preparation(qes::Program<>& program) {
    std::set<sptr<raw_vertex_t>> checks_this_stage = get_checks_at_stage(PREP_STAGE);
    // First, perform the H gates.
    std::vector<uint64_t> h_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        push_back_range(h_operands, prep_tear_get_h_operands(rpq, PREP_STAGE));
    }
    safe_emplace_back(program, "h", h_operands);
    // Now, we must perform the CNOTS.
    while (1) {
        cx_in_use_set.clear();
        std::vector<uint64_t> cx_operands;
        for (sptr<raw_vertex_t> rpq : checks_this_stage) {
            push_back_range(cx_operands, prep_tear_get_cx_operands(rpq, PREP_STAGE));
        }
        if (cx_operands.empty()) break;
        safe_emplace_back(program, "cx", cx_operands);
        // Perform any proxy resets:
        perform_proxy_resets(program);
    }
    // Now, we are done with preparation. We set stage_map to 1 if:
    //  (1) h_gate_stage_map of Z flags/X checks are set to PREP_STAGE, and
    //  (2) flag_stage_map of flags for each check is set to PREP_STAGE
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        auto support = get_support(rpq);
        // Check condition 1.
        bool cond_met = true;
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            cond_met &= (h_gate_stage_map[rpq] == PREP_STAGE);
        } else {
            for (sptr<raw_vertex_t> rfq : support.flags) {
                cond_met &= (h_gate_stage_map[rfq] == PREP_STAGE);
            }
        }
        for (sptr<raw_vertex_t> rfq : support.flags) {
            cond_met &= (flag_stage_map[rfq] == PREP_STAGE);
        }
        if (cond_met) stage_map[rpq] = BODY_STAGE;
    }
}

void
Scheduler::build_body(qes::Program<>& program) {
    data_cx_map.clear();

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
        check_operator_max_weight = std::max(get_support(rpq).data.size(), check_operator_max_weight);
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
            lp_var_t x = LP.add_var(rdq, 1, upper_bound, lp_var_t::bounds::both, lp_var_t::domain::integer);

            lp_constr_t max_con(max_of_all, x, lp_constr_t::direction::ge);
            LP.add_constraint(max_con);
        }
        // Add uniqueness constraints
        for (sptr<raw_vertex_t> rq1 : partial_support) {
            lp_var_t x = LP.get_var(rq1);
            for (sptr<raw_vertex_t> rq2 : partial_support) {
                if (rq1 >= rq2) continue;
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
                lp_constr_t con3(t - x, M*(1-y), lp_constr_t::direction::le);
                LP.add_constraint(con2);
                LP.add_constraint(con3);
                ind_sum += y;
                running_ind_ctr_map[rpq][rx]++;
                running_ind_ctr_map[rx][rpq]++;
            }
        }
        for (auto p : ind_sum_map) {
            // Check if the total number of crossing qubits is even. If so, add this constraint.
            if (running_ind_ctr_map[rpq][p.first] % 2 == 1) continue;

            lp_expr_t ind_sum = p.second + running_ind_sum_map[rpq][p.first];
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
            std::cerr << "when scheduling for " << print_v(rpq) << ":" << std::endl;
            for (auto rdq : partial_support) {
                std::cerr << "\t" << print_v(rdq) << ", t =";
                for (auto p : existing_data_cnot_times[rdq]) {
                    std::cerr << " " << p.first << "(" << print_v(p.second) << ")";
                }
                std::cerr << std::endl;
            }
            std::cerr << "possible problems:" << std::endl
                    << "\tnumber of variables (" << LP.variables.size() << ")" << std::endl
                    << "\tnumber of constraints (" << LP.l_constraints.size() << ")" << std::endl;
            exit(1);
        }
        // Update the schedule.
        std::vector<sptr<raw_vertex_t>> check_sch(static_cast<size_t>(round(obj)), nullptr);
        for (sptr<raw_vertex_t> rdq : partial_support) {
            size_t t = static_cast<size_t>(round(LP.get_value(rdq)));
            existing_data_cnot_times[rdq].push_back(std::make_pair(t, rpq));
            check_sch[t-1] = rdq;
        }
        data_cnot_schedules[rpq] = check_sch;
        // Update running_ind_sum_map.
        for (auto p : ind_sum_map) {
            sptr<raw_vertex_t> rx = p.first;
            lp_expr_t ind_sum = p.second;
            // Compute ind_sum value.
            double v = ind_sum.constant;
            for (auto pp : ind_sum.l_coefs) {
                v += pp.second * LP.get_value(pp.first);
            }
            running_ind_sum_map[rpq][rx] += v;
            running_ind_sum_map[rx][rpq] += v;
        }
    }
    // Now, that we have the CNOT schedules of all the partial supports, we need to execute all possible
    // CNOTs for this cycle.
    size_t k = 0;
    while (1) {
        cx_in_use_set.clear();
        std::vector<uint64_t> cx_operands;
        for (sptr<raw_vertex_t> rpq : checks_this_stage) {
            push_back_range(cx_operands, body_get_cx_operands(rpq, data_cnot_schedules[rpq], k));
        }
        if (cx_operands.empty()) break;
        safe_emplace_back(program, "cx", cx_operands);
        perform_proxy_resets(program);
        k++;
    }
    // We are done with the CX body. We advance a parity qubit if:
    //  (1) all data qubits have data_stage_map = BODY_STAGE.
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        auto& support = get_support(rpq);
        bool cond_met = true;
        for (sptr<raw_vertex_t> rdq : support.data) {
            cond_met &= (data_stage_map[rpq][rdq] == BODY_STAGE);
        }
        if (cond_met) stage_map[rpq] = TEAR_STAGE;
    }
}

void
Scheduler::build_teardown(qes::Program<>& program) {
    std::set<sptr<raw_vertex_t>> checks_this_stage = get_checks_at_stage(TEAR_STAGE);
    // As teardown is just the inversion of preparation, perfrom the CX gates first.
    while (1) {
        cx_in_use_set.clear();
        std::vector<uint64_t> cx_operands;
        for (sptr<raw_vertex_t> rpq : checks_this_stage) {
            push_back_range(cx_operands, prep_tear_get_cx_operands(rpq, TEAR_STAGE));
        }
        if (cx_operands.empty()) break;
        safe_emplace_back(program, "cx", cx_operands);
        perform_proxy_resets(program);
    }
    // Perform the H gates.
    //
    // However, note that since the H gates are dependent on the CX gates, we need to
    // check if the CX gates are complete for the flags and parity qubits.
    std::vector<uint64_t> h_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        bool all_cx_done = true;
        for (sptr<raw_vertex_t> rfq : get_support(rpq).flags) {
            all_cx_done &= (flag_stage_map[rfq] >= TEAR_STAGE);
        }
        if (all_cx_done) {
            push_back_range(h_operands, prep_tear_get_h_operands(rpq, TEAR_STAGE));
        }
    }
    safe_emplace_back(program, "h", h_operands);
    // Now, we are done with teardown. We advance the stage if
    //  (1) h_gate_state_map = TEAR_STAGE for the X-check/Z-flags
    //  (2) flag_stage_map = TEAR_STAGE.
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        auto support = get_support(rpq);
        // Check condition 1.
        bool cond_met = true;
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            cond_met &= (h_gate_stage_map[rpq] == TEAR_STAGE);
        } else {
            for (sptr<raw_vertex_t> rfq : support.flags) {
                cond_met &= (h_gate_stage_map[rfq] == TEAR_STAGE);
            }
        }
        for (sptr<raw_vertex_t> rfq : support.flags) {
            cond_met &= (flag_stage_map[rfq] == TEAR_STAGE);
        }
        if (cond_met) stage_map[rpq] = MEAS_STAGE;
    }
}

void
Scheduler::build_measurement(qes::Program<>& program) {
    std::vector<uint64_t> m_operands;

    std::set<sptr<raw_vertex_t>> checks_this_stage = get_checks_at_stage(MEAS_STAGE);
    // Measure all parity qubits and flag qubits.
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        push_back_measurement(m_operands, rpq);
        for (sptr<raw_vertex_t> rfq : get_support(rpq).flags) {
            push_back_measurement(m_operands, rfq);
        }
    }
    safe_emplace_back(program, "measure", m_operands);
    // Now, we need to reset all parity, flags, and proxies.
    // We also need to clear out the corresponding active_role_map entries.
    std::vector<uint64_t> r_operands(m_operands);
    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        auto& support = get_support(rpq);
        // Reset the proxies.
        // 
        // Also, we should release all parity qubits, flag qubits, and proxy qubits.
        for (sptr<raw_vertex_t> rprx : support.proxies) {
            sptr<phys_vertex_t> pprx = net_p->role_to_phys[rprx];
            r_operands.push_back(pprx->id);
            // Release the physical qubit for the proxy.
            release_physical_qubit(rprx);
        }
        // Now do flags:
        for (sptr<raw_vertex_t> rfq : support.flags) {
            release_physical_qubit(rfq);
        }
        release_physical_qubit(rpq);
#ifdef PROTEAN_DEBUG
        std::cerr << "[ build_measurement ] " << print_v(rpq) << " has ended @ cycle = " << cycle << "\n";
#endif
    }
    safe_emplace_back(program, "reset", r_operands);

    for (sptr<raw_vertex_t> rpq : checks_this_stage) {
        stage_map[rpq] = DONE_STAGE;
    }
}

std::vector<uint64_t>
Scheduler::prep_tear_get_h_operands(sptr<raw_vertex_t> rpq, stage_t s) {
    std::vector<uint64_t> operands;
    std::vector<sptr<raw_vertex_t>> _operands;

    auto& support = get_support(rpq);
    // Collect potential operands.
    if (rpq->qubit_type == raw_vertex_t::type::xparity) {
        _operands.push_back(rpq);
    } else {
        push_back_range(_operands, support.flags);
    }
    // Perform the H gates if possible.
    for (sptr<raw_vertex_t> rx : _operands) {
        // Check if we have already done this gate in this stage.
        if (h_gate_stage_map[rx] >= s) continue;
        sptr<phys_vertex_t> px = net_p->role_to_phys[rx];
        if (is_good_for_current_cycle(rx)) {
            operands.push_back(px->id);
            h_gate_stage_map[rx] = s;
        }
    }
    return operands;
}

std::vector<uint64_t>
Scheduler::prep_tear_get_cx_operands(sptr<raw_vertex_t> rpq, Scheduler::stage_t s) {
    std::vector<uint64_t> operands;
    auto& support = get_support(rpq);
    for (sptr<raw_vertex_t> rfq : support.flags) {
        // Check if the flag is already done.
        if (flag_stage_map[rfq] >= s) continue;
        // Attempt to do the flag.
        cx_t cx = get_next_edge_between(rfq, rpq, rpq->qubit_type == raw_vertex_t::type::xparity, s);
        // Check CX status.
        if (cx_return_status) {
            // Something went wrong and we could not get an edge.
            //
            // Check if the edge is all done:
            if (cx_return_status == CX_RET_FINISHED) flag_stage_map[rfq] = s;
            continue;
        }
        push_back_cx(operands, cx, s);
    }
    return operands;
}

std::set<sptr<raw_vertex_t>>
Scheduler::body_get_partial_data_support(sptr<raw_vertex_t> rpq) {
    const bool is_x_check = (rpq->qubit_type == raw_vertex_t::type::xparity);
    std::set<sptr<raw_vertex_t>> partial_support;

    RawNetwork& raw_net = net_p->raw_connection_network;

    auto& support = get_support(rpq);
    for (sptr<raw_vertex_t> rdq : support.data) {
        // Four possibilites:
        //  (1) rdq is directly connected to rpq (no flag)
        //  (2) rdq is indirectly connected to rpq (proxy walk)
        //  (3) rdq is directly connected to a flag
        //  (4) rdq is indirectly connected to a flag
        sptr<raw_vertex_t> other = nullptr;
        if (raw_net.contains(rdq, rpq)) {
            // Situation 1
            other = rpq;
        } else if (raw_net.flag_assignment_map[rpq].count(rdq)) {
            sptr<raw_vertex_t> rfq = raw_net.flag_assignment_map[rpq][rdq];
            // Check if rdq is connected to rfq.
            if (raw_net.contains(rdq, rfq)) {
                // Situation 3
                other = rfq;
            } else {
                // Situation 4.
                auto& proxy_walk_path = is_x_check 
                                        ? get_proxy_walk_path(rfq, rdq)
                                        : get_proxy_walk_path(rdq, rfq);
                other = test_and_get_other_endpoint_if_ready(rdq, proxy_walk_path);
            }
        } else {
            // Situation 2.
            auto& proxy_walk_path = is_x_check
                                    ? get_proxy_walk_path(rpq, rdq) 
                                    : get_proxy_walk_path(rdq, rpq);
            other = test_and_get_other_endpoint_if_ready(rdq, proxy_walk_path);
        }
        // If other == nullptr, this means that we are not ready to do the edge in the proxy walk.
        if (other == nullptr) continue;
        if (!is_good_for_current_cycle(rdq) || !is_good_for_current_cycle(other)) continue;
        sptr<raw_edge_t> re = raw_net.get_edge(rdq, other);
        if (visited_edge_map[re] < BODY_STAGE) {
            partial_support.insert(rdq);
            // Record the CX gate for future use.
            cx_t cx;
            if (is_x_check) {
                cx = std::make_tuple(other, rdq, re);
            } else {
                cx = std::make_tuple(rdq, other, re);
            }
            tlm_put(data_cx_map, rpq, rdq, cx);
        }
    }
    return partial_support;
}

std::vector<uint64_t>
Scheduler::body_get_cx_operands(
        sptr<raw_vertex_t> rpq,
        const std::vector<sptr<raw_vertex_t>>& cnot_sch,
        size_t k)
{
    RawNetwork& raw_net = net_p->raw_connection_network;

    std::vector<uint64_t> operands;

    // Check if there is a data qubit at time step k.
    sptr<raw_vertex_t> s_rdq = k < cnot_sch.size() ? cnot_sch.at(k) : nullptr;
    // We are obligated to first perform this CX.
    if (s_rdq != nullptr) {
        cx_t cx = data_cx_map[rpq][s_rdq];
        push_back_cx(operands, cx, BODY_STAGE);
    }
    // Now, we can look to performing other data qubit CNOTs.
    for (sptr<raw_vertex_t> rdq : get_support(rpq).data) {
        // Check if the data qubit is already done.
        if (data_stage_map[rpq][rdq] >= BODY_STAGE || s_rdq == rdq) continue;
        // Check if rdq is connected to a flag or rpq.
        sptr<raw_vertex_t> other = rpq;
        if (raw_net.flag_assignment_map[rpq].count(rdq)) {
            other = raw_net.flag_assignment_map[rpq][rdq];
        }
        cx_t cx = get_next_edge_between(rdq, other, rpq->qubit_type == raw_vertex_t::type::xparity, BODY_STAGE);
        if (cx_return_status) {
            // Something went wrong and we could not get an edge.
            //
            // Check if the edge is all done:
            if (cx_return_status == CX_RET_FINISHED) {
#ifdef PROTEAN_DEBUG
                std::cerr << "finished CX(" << print_v(rdq) << ", " << print_v(other) << ")" << std::endl;
#endif
                data_stage_map[rpq][rdq] = BODY_STAGE;
            }
            continue;
        }
        // Check if one of the endpoints of the cx is a data qubit. If so, do not do the CX -- it needs
        // to wait its turn.
        sptr<raw_vertex_t> rx = std::get<0>(cx),
                            ry = std::get<1>(cx);
        if (rx == rdq || ry == rdq) continue;
        push_back_cx(operands, cx, BODY_STAGE);
    }
    return operands;
}

qes::Program<>
PhysicalNetwork::make_schedule() {
    qes::Program<> program;
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
    Scheduler scheduler(this);
    // Set up h_operands (data qubits, used if is_memory_x), r_operands (all qubits), m_operands (data qubits).
    std::map<sptr<raw_vertex_t>, size_t> data_meas_ctr_map;
    std::vector<uint64_t> h_operands, r_operands, m_operands;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        r_operands.push_back(pv->id);
    }
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        if (rv->qubit_type != raw_vertex_t::type::data) continue;
        sptr<phys_vertex_t> pv = role_to_phys[rv];
        data_meas_ctr_map[rv] = m_operands.size();
        m_operands.push_back(pv->id);
        if (config.is_memory_x) {
            h_operands.push_back(pv->id);
        }
    }
    // Prologue:
    safe_emplace_back(program, "reset", r_operands);
    safe_emplace_back(program, "h", h_operands);
    // Syndrome Extraction:
    qes::Program<> round = scheduler.run();

    const size_t n_mt = scheduler.get_measurement_ctr();
    size_t event_ctr = 0;
    size_t meas_ctr_offset = 0;
    for (size_t r = 0; r < config.rounds; r++) {
        push_back_range(program, round);
        // Build detection events.
        for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
            bool is_x_check = (rv->qubit_type == raw_vertex_t::type::xparity);
            bool is_z_check = (rv->qubit_type == raw_vertex_t::type::zparity);
            if (!is_z_check && !is_x_check) continue;

            if (is_x_check == config.is_memory_x) {
                const size_t mt = scheduler.get_measurement_time(rv) + meas_ctr_offset;
                // Make detection event for the parity qubit.
                std::vector<uint64_t> operands{event_ctr, mt};
                if (r > 0) operands.push_back(mt - n_mt);
                safe_emplace_back(program, "event", operands);
//              program.back().put(print_v(rv));
                event_ctr++;
            } else {
                // Make detection events for the flag qubits.
                for (sptr<raw_vertex_t> rfq : raw_connection_network.flag_ownership_map[rv]) {
                    const size_t mt = scheduler.get_measurement_time(rfq) + meas_ctr_offset;
                    safe_emplace_back(program, "event", {event_ctr, mt});
//                  program.back().put(print_v(rfq));
                    event_ctr++;
                }
            }
        }
        meas_ctr_offset += n_mt;
    }
    // Epilogue:
    const size_t prev_meas_ctr_offset = meas_ctr_offset - n_mt;
    safe_emplace_back(program, "h", h_operands);
    safe_emplace_back(program, "measure", m_operands);
    // Build detection events and observable.
    TannerGraph& tanner_graph = raw_connection_network.tanner_graph;
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        bool is_x_check = (rv->qubit_type == raw_vertex_t::type::xparity);
        bool is_z_check = (rv->qubit_type == raw_vertex_t::type::zparity);
        if (!is_z_check && !is_x_check) continue;
        if (is_x_check != config.is_memory_x) continue;

        const size_t pq_mt = scheduler.get_measurement_time(rv) + prev_meas_ctr_offset;
        std::vector<uint64_t> operands{event_ctr, pq_mt};

        sptr<tanner::vertex_t> tv = raw_connection_network.v_tanner_raw_map.at(rv);
        for (sptr<tanner::vertex_t> tdq : tanner_graph.get_neighbors(tv)) {
            sptr<raw_vertex_t> rdq = raw_connection_network.v_tanner_raw_map.at(tdq);
            operands.push_back(data_meas_ctr_map[rdq] + meas_ctr_offset);
        }
        safe_emplace_back(program, "event", operands);
//      program.back().put(print_v(rv));
        event_ctr++;
    }
    auto obs_list = tanner_graph.get_obs(config.is_memory_x);
    size_t obs_ctr = 0;
    for (auto& obs : obs_list) {
        std::vector<uint64_t> operands{obs_ctr};
        for (sptr<tanner::vertex_t> tdq : obs) {
            sptr<raw_vertex_t> rdq = raw_connection_network.v_tanner_raw_map.at(tdq);
            operands.push_back(data_meas_ctr_map[rdq] + meas_ctr_offset);
        }
        safe_emplace_back(program, "obs", operands);
        obs_ctr++;
    }
    return program;
}

sptr<raw_vertex_t>
Scheduler::test_and_get_other_endpoint_if_ready(
        sptr<raw_vertex_t> endpoint,
        std::vector<sptr<raw_vertex_t>> path) 
{
    /*
    auto src = path[0],
            dst = path[path.size()-1];
    for (size_t i = 1; i < path.size(); i++) {
        sptr<raw_vertex_t> rx = path[i-1],
                            ry = path[i];
        sptr<raw_edge_t> re = net_p->raw_connection_network.get_edge(rx, ry);
        if (rx == endpoint || ry == endpoint) {
            if (visited_edge_map[re] < BODY_STAGE) {
                // Check if the edge is even doable.
                if (has_contention(rx) || has_contention(ry)) return nullptr;
                if (!is_good_for_current_cycle(rx) || !is_good_for_current_cycle(ry)) return nullptr;
                if (!is_proxy_usable(rx, src, dst) || !is_proxy_usable(ry, src, dst)) return nullptr;
                // Then the edge is not done.
                return rx == endpoint ? ry : rx;
            }
        } else {
            if (visited_edge_map[re] < BODY_STAGE) {
                return nullptr;
            }
        }
    }
    */
    cx_t cx = get_next_edge_between(path[0], path.back(), false, BODY_STAGE);
    if (cx_return_status) return nullptr;

    sptr<raw_vertex_t> rx = std::get<0>(cx),
                        ry = std::get<1>(cx);
    if (rx == endpoint)         return ry;
    else if (ry == endpoint)    return rx;
    else                        return nullptr;
}

Scheduler::cx_t
Scheduler::get_next_edge_between(sptr<raw_vertex_t> src, sptr<raw_vertex_t> dst, bool for_x_check, stage_t s) {
    if (for_x_check) return get_next_edge_between(dst, src, false, s);

    cx_return_status = CX_RET_GOOD;
    // Check if src and dst are directly connected, or are connected via proxy.
    RawNetwork& raw_net = net_p->raw_connection_network;
    if (raw_net.contains(src, dst)) {
        sptr<raw_edge_t> re = raw_net.get_edge(src, dst);
        if (visited_edge_map[re] < s) {
            if (has_contention(src) || has_contention(dst)) {
                return ret_null_and_set_status(CX_RET_CONTENTION);
            }
            if (!is_good_for_current_cycle(src) || !is_good_for_current_cycle(dst)) {
                return ret_null_and_set_status(CX_RET_TOO_EARLY);
            }
            return std::make_tuple(src, dst, re);
        }
    } else {
        // Otherwise, they are connected via proxy.
        std::vector<sptr<raw_vertex_t>>& path = get_proxy_walk_path(src, dst);
        for (size_t i = 1; i < path.size(); i++) {
            sptr<raw_vertex_t> rx = path[i-1],
                                ry = path[i];
            sptr<raw_edge_t> re = raw_net.get_edge(rx, ry);
            if (visited_edge_map[re] < s) {
                if (has_contention(rx) || has_contention(ry)) {
#ifdef PROTEAN_DEBUG
                    std::cerr << "CX(" << print_v(src) << ", " << print_v(dst) << ") has contention: (" 
                        << print_v(path[i-1]) << ", " << print_v(path[i]) << ") @ cycle = " << cycle << std::endl;
#endif
                    return ret_null_and_set_status(CX_RET_CONTENTION);
                }
                if (!is_good_for_current_cycle(rx) || !is_good_for_current_cycle(ry)) {
#ifdef PROTEAN_DEBUG
                    sptr<phys_vertex_t> pv = net_p->role_to_phys[path[i-1]],
                                        pw = net_p->role_to_phys[path[i]];
                    std::cerr << "CX(" << print_v(src) << ", " << print_v(dst) << ") cannot proceed: " 
                        << "need " << print_v(path[i-1]) << " to reach " << pv->cycle_role_map.at(path[i-1])
                        << " and " << print_v(path[i]) << " to reach " << pw->cycle_role_map.at(path[i])
                        << ", or need " << print_v(active_role_map[pv])
                        << " and/or " << print_v(active_role_map[pw])
                        << " to release resource @ cycle = " << cycle << "\n";
#endif
                    return ret_null_and_set_status(CX_RET_TOO_EARLY);
                }
                if (!is_proxy_usable(rx, src, dst) || !is_proxy_usable(ry, src, dst)) {
#ifdef PROTEAN_DEBUG
                    auto& erx = proxy_occupied_map[rx];
                    auto& ery = proxy_occupied_map[ry];
                    std::cerr << print_v(src) << " and " << print_v(dst) << " cannot acquire proxies: "
                        << print_v(rx) << " [ " << print_v(std::get<0>(erx)) << ", "
                        << print_v(std::get<1>(erx)) << ", " << std::get<2>(erx) << " ] and "
                        << print_v(ry) << " [ " << print_v(std::get<0>(ery)) << ", "
                        << print_v(std::get<1>(ery)) << ", " << std::get<2>(ery) << " ]" << std::endl;
#endif
                    return ret_null_and_set_status(CX_RET_PROXY_OCCUPIED);
                }
                // Check if the first operand in the path is a proxy. If so, then mark for reset (it is done
                // transfering state).
                if (rx->qubit_type == raw_vertex_t::type::proxy) {
                    proxy_reset_map[rx] = 1;
#ifdef PROTEAN_DEBUG
                    std::cerr << print_v(src) << " and " << print_v(dst) << " will reset " << print_v(rx)
                        << " for CX(" << print_v(rx) << ", " << print_v(ry) << ")" << std::endl;
#endif
                }
                // Similarly, set proxy_occupied_map as well.
                if (ry->qubit_type == raw_vertex_t::type::proxy) {
                    proxy_occupied_map[ry] = std::make_tuple(src, dst, 1);
#ifdef PROTEAN_DEBUG
                    std::cerr << print_v(src) << " and " << print_v(dst) << " acquired " << print_v(ry)
                        << " for CX(" << print_v(rx) << ", " << print_v(ry) << ")" << std::endl;
#endif
                }
                // Then, return this edge -- it is not done.
                return std::make_tuple(rx, ry, re);
            }
        }
    }
    // No edge to be done:
    return ret_null_and_set_status(CX_RET_FINISHED);
}

}   // protean
}   // qontra
