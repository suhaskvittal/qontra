/*
 *  author: Suhas Vittal
 *  date:   7 January 2024
 * */

#include "qontra/protean/scheduler.h"

#include <vtils/linprog/manager.h>
#include <vtils/two_level_map.h>
#include <vtils/utility.h>

#include <deque>
#include <iostream>
#include <sstream>

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;
using namespace vtils;

template <class PTR> std::string
print_path_inl(std::vector<PTR> path) {
    std::ostringstream sout;
    for (auto v : path) {
        sout << " " << print_v(v);
    }
    return sout.str();
}

inline void
safe_emplace_back(qes::Program<>& program, std::string name, std::vector<int64_t> operands) {
    if (operands.empty()) return;
    program.emplace_back(name, std::vector<int64_t>(operands.begin(), operands.end()));
}

inline void
safe_emplace_back(qes::Program<>& program, std::string name, std::vector<uint64_t> operands) {
    if (operands.empty()) return;
    program.emplace_back(name, std::vector<int64_t>(operands.begin(), operands.end()));
}

Scheduler::Scheduler(PhysicalNetwork* n_p) 
    :active_role_map(),
    checks_this_cycle(),
    meas_ctr_map(),
    cycle(0),
    mctr(0),
    tanner_graph(n_p->get_tanner_graph()),
    raw_network(n_p->get_raw_connection_network()),
    network(n_p)
{}

qes::Program<>
Scheduler::run() {
    qes::Program<> program;
    std::set<sptr<raw_vertex_t>> finished_checks;
    cycle = 0;
    while (true) {
        std::cout << "====== CYCLE " << cycle << " ======" << std::endl;
        std::cout << "Checks this cycle:";

        size_t k = program.size();

        checks_this_cycle.clear();
        for (sptr<raw_vertex_t> rv : raw_network->get_vertices()) {
            if (rv->qubit_type == raw_vertex_t::type::xparity
                || rv->qubit_type == raw_vertex_t::type::zparity)
            {
                sptr<phys_vertex_t> pv = network->role_to_phys.at(rv);
                size_t cyc = pv->cycle_role_map.at(rv);
                if (cyc <= cycle && !finished_checks.count(rv)) {
                    checks_this_cycle.push_back(rv);

                    std::cout << " " << print_v(rv);
                }
            }
        }
        std::cout << std::endl;

        if (checks_this_cycle.empty()) break;
        build_preparation(program);
        build_body(program);
        build_teardown(program);
        
        for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
            finished_checks.insert(rpq);
        }
        cycle++;
        // Inject timing error on first instruction of this cycle.
        if (k < program.size()) {
            program[k].put("timing_error");
        }
    }
    return program;
}

void
Scheduler::build_preparation(qes::Program<>& program) {
    prep_tear_h_gates(program);
    prep_tear_cx_gates(program);
}

void
Scheduler::build_body(qes::Program<>& program) {
    auto cx_sch_map = compute_schedules();
    // Perform the CX gates. Note that the schedule contains times
    // for data CX to their target (either parity or flag).
    size_t k = 0;
    bool any_cx_done;
    do {
        std::vector<cx_t> cx_arr;
        any_cx_done = false;
        for (const auto& p : cx_sch_map) {
            sptr<raw_vertex_t> rsq = p.first;
            const std::vector<sptr<raw_vertex_t>>& sch = p.second;
            // Make sure rsq is valid.
            test_and_set_exit_on_fail(rsq, "build_body");
            // Now get the k-th data qubit.
            sptr<raw_vertex_t> rdq = sch.at(k);
            if (rdq == nullptr) continue;
            sptr<raw_vertex_t> common_rpq = raw_network->are_in_same_support(rdq, rsq);
            if (common_rpq->qubit_type == raw_vertex_t::type::xparity) {
                cx_arr.emplace_back(rsq, rdq);
            } else {
                cx_arr.emplace_back(rdq, rsq);
            }
            any_cx_done = true;
        }
        schedule_cx_along_path(cx_arr, program);
        k++;
    } while (any_cx_done);
}

void
Scheduler::build_teardown(qes::Program<>& program) {
    prep_tear_cx_gates(program);
    prep_tear_h_gates(program);
    // Measure the parity and flag qubits. Reset all proxies.
    std::vector<int64_t> m_operands, r_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        push_back_measurement(m_operands, rpq);
        r_operands.push_back(qu(rpq));
        release_qubit(rpq);

        auto& support = get_support(rpq);
        for (sptr<raw_vertex_t> rfq : support.flags) {
            push_back_measurement(m_operands, rfq);
            r_operands.push_back(qu(rfq));
            release_qubit(rfq);
        }
        for (sptr<raw_vertex_t> rprx : support.proxies) {
            r_operands.push_back(qu(rprx));
        }
    }
    safe_emplace_back(program, "measure", m_operands);
    safe_emplace_back(program, "reset", r_operands);
}

void
Scheduler::prep_tear_h_gates(qes::Program<>& program) {
    std::vector<int64_t> h_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        test_and_set_exit_on_fail(rpq, "build_preparation");
        // Perform an H gate if necessary.
        auto& support = get_support(rpq);
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            h_operands.push_back(qu(rpq));
        } else {
            for (sptr<raw_vertex_t> rfq : support.flags) {
                test_and_set_exit_on_fail(rfq, "build_preparation");
                h_operands.push_back(qu(rfq));
            }
        }
    }
    safe_emplace_back(program, "h", h_operands);
}

void
Scheduler::prep_tear_cx_gates(qes::Program<>& program) {
    std::vector<cx_t> cx_arr;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        auto& support = get_support(rpq);
        for (sptr<raw_vertex_t> rfq : support.flags) {
            if (rpq->qubit_type == raw_vertex_t::type::xparity) {
                cx_arr.emplace_back(rpq, rfq);
            } else {
                cx_arr.emplace_back(rfq, rpq);
            }
        }
    }
    schedule_cx_along_path(cx_arr, program);
}

void
Scheduler::schedule_cx_along_path(const std::vector<cx_t>& cx_arr, qes::Program<>& program) {
    std::vector<std::vector<sptr<raw_vertex_t>>> path_arr(cx_arr.size());
    std::vector<size_t> k_arr(cx_arr.size(), 1);

    for (size_t i = 0; i < cx_arr.size(); i++) {
        cx_t cx = cx_arr.at(i);
        sptr<net::raw_vertex_t> src = std::get<0>(cx),
                                dst = std::get<1>(cx);
        path_arr[i] = raw_network->get_proxy_walk_path(src, dst);
    }
    // Also track undo paths:
    std::vector<std::vector<sptr<raw_vertex_t>>> undo_arr;
    std::vector<size_t> undo_k_arr;
    // Now perform CX along these paths.
    std::set<sptr<net::raw_vertex_t>> proxies_in_use;
    while (true) {
        std::set<int64_t> in_use;
        std::vector<int64_t> cx_operands;
        for (size_t i = 0; i < cx_arr.size(); i++) {
            const auto& path = path_arr[i];
            size_t& k = k_arr[i];
            if (!try_and_push_back_cx_operands(
                    cx_operands,
                    in_use,
                    path,
                    k,
                    [&] (sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry)
                    {
                        return !proxies_in_use.count(ry);
                    })
                )
            {
                continue;
            }

            sptr<raw_vertex_t> ry = path.at(k);
            if (ry->qubit_type == raw_vertex_t::type::proxy) {
                proxies_in_use.insert(ry);
            }
            // Check if we have a valid undo path. If so, mark it as a path
            // we need to complete.
            k++;
            if (k == path.size() && path.size() > 2) {
                undo_arr.emplace_back(path.rbegin()+1, path.rend());
                undo_k_arr.push_back(1);
            }
        }
        // Now do the undo paths.
        for (size_t i = 0; i < undo_arr.size(); i++) {
            const auto& path = undo_arr[i];
            size_t& k = undo_k_arr[i];
            if (!try_and_push_back_cx_operands(
                    cx_operands,
                    in_use,
                    path,
                    k,
                    [&] (sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry) { return true; },
                    true)
                )
            {
                continue;
            }

            sptr<raw_vertex_t> rx = path.at(k-1);
            if (rx->qubit_type == raw_vertex_t::type::proxy) {
                proxies_in_use.erase(rx);
            }
            k++;
        }
        if (cx_operands.empty()) break;
        safe_emplace_back(program, "cx", cx_operands);
    }
}

std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>>
Scheduler::compute_schedules() {
    std::map<sptr<raw_vertex_t>, std::vector<sptr<raw_vertex_t>>> sch_map;
    // Compute max operator weight.
    size_t max_operator_weight = 0;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        auto& support = get_support(rpq);
        max_operator_weight = std::max(max_operator_weight, support.data.size());
    }
    const fp_t upper_bound = static_cast<fp_t>(2 * max_operator_weight);
    // Now, compute schedule using LPs:
    std::map<sptr<raw_vertex_t>, std::vector<std::tuple<size_t, sptr<raw_vertex_t>>>>
        existing_cnot_times;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        auto& support = get_support(rpq);
        // Make the LP:
        CPXLPManager<sptr<raw_vertex_t>> LP;
        lp_var_t max_of_all = 
            LP.add_slack_var(0, upper_bound, lp_var_t::bounds::both, lp_var_t::domain::integer);
        for (sptr<raw_vertex_t> rdq : support.data) {
            lp_var_t x = LP.add_var(rdq, 1, upper_bound, lp_var_t::bounds::both, lp_var_t::domain::integer);
            lp_constr_t max_con(max_of_all, x, lp_constr_t::direction::ge);
            LP.add_constraint(max_con);
        }
        // Add uniqueness constraints.
        for (sptr<raw_vertex_t> rx : support.data) {
            lp_var_t x = LP.get_var(rx);
            sptr<raw_vertex_t> rx_rsq =
                raw_network->flag_assignment_map[rpq].count(rx) ? raw_network->flag_assignment_map[rpq][rx] : rpq;
            for (sptr<raw_vertex_t> ry : support.data) {
                if (rx >= ry) continue;
                lp_var_t y = LP.get_var(ry);
                sptr<raw_vertex_t> ry_rsq =
                    raw_network->flag_assignment_map[rpq].count(ry) ? raw_network->flag_assignment_map[rpq][ry] : rpq;
                // Only add the following constraint if rx and ry both interact
                // with rpq or the same flag.
                if (rx_rsq == ry_rsq) {
                    lp_constr_t con(x, y, lp_constr_t::direction::neq);
                    LP.add_constraint(con);
                }
            }
        }
        // Add stabilizer constraints:
        std::map<sptr<raw_vertex_t>, lp_expr_t> ind_sum_map;
        for (sptr<raw_vertex_t> rdq : support.data) {
            lp_var_t x = LP.get_var(rdq);
            for (auto& tt : existing_cnot_times[rdq]) {
                size_t t = std::get<0>(tt);
                sptr<raw_vertex_t> owner = std::get<1>(tt);
                // Add conflict constraint.
                lp_constr_t con1(x, t, lp_constr_t::direction::neq);
                LP.add_constraint(con1);
                // if is_x_check != tt_from_x, then we need to add commutation
                // constraints.
                if (rpq->qubit_type == owner->qubit_type) continue;
                lp_expr_t& ind_sum = ind_sum_map[owner];

                const double M = 100'000;
                lp_var_t y = LP.add_slack_var(0, 1, lp_var_t::bounds::both, lp_var_t::domain::binary);
                // x - t <= My and t - x <= M(1-y).
                lp_constr_t con2(x - t, M*y, lp_constr_t::direction::le),
                            con3(t - x, M*(1-y), lp_constr_t::direction::le);
                LP.add_constraint(con2);
                LP.add_constraint(con3);

                ind_sum += y;
            }
        }
        for (auto& p : ind_sum_map) {
            // Constrain ind_sum to be even.
            lp_expr_t& ind_sum = p.second;
            lp_var_t y = LP.add_slack_var(0, 0, lp_var_t::bounds::lower, lp_var_t::domain::integer);
            lp_constr_t con(ind_sum, 2*y, lp_constr_t::direction::eq);
            LP.add_constraint(con);
        }
        // Now, solve the LP.
        LP.build(max_of_all, false);
        double obj;
        int solstat, status;
        if ((status=LP.solve(&obj, &solstat))) {
            std::cerr << "build_body: program is infeasible, solstat = " << solstat
                << ", status = " << status << std::endl;
            std::cerr << "when scheduling for " << print_v(rpq) << ":" << std::endl;
            for (auto rdq : support.data) {
                std::cerr << "\t" << print_v(rdq) << ", t =";
                for (auto tt : existing_cnot_times[rdq]) {
                    std::cerr << " " << std::get<0>(tt) << "(" << print_v(std::get<1>(tt)) << ")";
                }
                std::cerr << std::endl;
            }
            std::cerr << "possible problems:" << std::endl
                    << "\tnumber of variables (" << LP.variables.size() << ")" << std::endl
                    << "\tnumber of constraints (" << LP.l_constraints.size() << ")" << std::endl;
            exit(1);
        }
        for (sptr<raw_vertex_t> rdq : support.data) {
            sptr<raw_vertex_t> rsq = 
                raw_network->flag_assignment_map[rpq].count(rdq) ? raw_network->flag_assignment_map[rpq][rdq] : rpq;
            size_t t = static_cast<size_t>(round(LP.get_value(rdq)));
            existing_cnot_times[rdq].emplace_back(t, rpq);
            if (!sch_map.count(rsq)) {
                sch_map[rsq] = std::vector<sptr<raw_vertex_t>>(upper_bound, nullptr);
            }
            sch_map[rsq][t-1] = rdq;
        }
    }
    return sch_map;
}

qes::Program<>
PhysicalNetwork::make_schedule() {
    qes::Program<> program;

    Scheduler scheduler(this);
    // Set up h_operands (data qubits, used if is_memory_x), r_operands (all qubits), m_operands (data qubits).
    std::map<sptr<raw_vertex_t>, size_t> data_meas_ctr_map;
    std::vector<uint64_t> h_operands, r_operands, m_operands;
    for (sptr<phys_vertex_t> pv : get_vertices()) {
        r_operands.push_back(pv->id);
    }
    for (sptr<raw_vertex_t> rv : raw_connection_network->get_vertices()) {
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
        for (sptr<raw_vertex_t> rv : raw_connection_network->get_vertices()) {
            bool is_x_check = (rv->qubit_type == raw_vertex_t::type::xparity);
            bool is_z_check = (rv->qubit_type == raw_vertex_t::type::zparity);
            if (!is_z_check && !is_x_check) continue;

            if (is_x_check == config.is_memory_x) {
                const size_t mt = scheduler.get_measurement_time(rv) + meas_ctr_offset;
                // Make detection event for the parity qubit.
                std::vector<uint64_t> operands{event_ctr, mt};
                if (r > 0) operands.push_back(mt - n_mt);
                safe_emplace_back(program, "event", operands);
                program.back().put(print_v(rv));
                event_ctr++;
            } else {
                // Make detection events for the flag qubits.
                for (sptr<raw_vertex_t> rfq : raw_connection_network->flag_ownership_map[rv]) {
                    const size_t mt = scheduler.get_measurement_time(rfq) + meas_ctr_offset;
                    safe_emplace_back(program, "event", std::vector<uint64_t>{event_ctr, mt});
                    program.back().put(print_v(rfq));
                    program.back().put("flag");
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
    for (sptr<raw_vertex_t> rv : raw_connection_network->get_vertices()) {
        bool is_x_check = (rv->qubit_type == raw_vertex_t::type::xparity);
        bool is_z_check = (rv->qubit_type == raw_vertex_t::type::zparity);
        if (!is_z_check && !is_x_check) continue;
        if (is_x_check != config.is_memory_x) continue;

        const size_t pq_mt = scheduler.get_measurement_time(rv) + prev_meas_ctr_offset;
        std::vector<uint64_t> operands{event_ctr, pq_mt};

        sptr<tanner::vertex_t> tv = raw_connection_network->v_tanner_raw_map.at(rv);
        for (sptr<tanner::vertex_t> tdq : tanner_graph->get_neighbors(tv)) {
            sptr<raw_vertex_t> rdq = raw_connection_network->v_tanner_raw_map.at(tdq);
            operands.push_back(data_meas_ctr_map[rdq] + meas_ctr_offset);
        }
        safe_emplace_back(program, "event", operands);
        program.back().put(print_v(rv));
        event_ctr++;
    }
    auto obs_list = tanner_graph->get_obs(config.is_memory_x);
    size_t obs_ctr = 0;
    for (auto& obs : obs_list) {
        std::vector<uint64_t> operands{obs_ctr};
        for (sptr<tanner::vertex_t> tdq : obs) {
            sptr<raw_vertex_t> rdq = raw_connection_network->v_tanner_raw_map.at(tdq);
            operands.push_back(data_meas_ctr_map[rdq] + meas_ctr_offset);
        }
        safe_emplace_back(program, "obs", operands);
        obs_ctr++;
    }
    return program;
}

}   // protean
}   // qontra
