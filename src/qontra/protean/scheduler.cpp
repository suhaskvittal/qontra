/*
 *  author: Suhas Vittal
 *  date:   7 January 2024
 * */

#include "qontra/protean/scheduler.h"

#include <vtils/linprog/manager.h>
#include <vtils/set_algebra.h>
#include <vtils/two_level_map.h>
#include <vtils/utility.h>

#include <deque>
#include <iostream>
#include <sstream>

#ifdef PROTEAN_PERF
#include <vtils/timer.h>
#endif

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
    program.emplace_back(name, operands);
}

inline void
safe_emplace_back(qes::Program<>& program, std::string name, std::set<int64_t> operands) {
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
{
    raw_network->enable_memoization = true;
}

qes::Program<>
Scheduler::run(size_t rounds, bool is_memory_x) {
    qes::Program<> program;
    //
    // PROLOGUE
    //
    std::vector<int64_t> r_operands, h_operands, m_operands;
    std::map<sptr<raw_vertex_t>, int64_t> data_meas_order_map;
    for (sptr<phys_vertex_t> pv : network->get_vertices()) {
        int64_t q = static_cast<int64_t>(pv->id);
        r_operands.push_back(q);
    }
    for (sptr<raw_vertex_t> rv : raw_network->get_vertices()) {
        if (rv->qubit_type != raw_vertex_t::type::data) continue;

        data_meas_order_map[rv] = h_operands.size();
        h_operands.push_back(qu(rv));
        m_operands.push_back(qu(rv));
    }
    safe_emplace_back(program, "reset", r_operands);
    if (is_memory_x) {
        safe_emplace_back(program, "h", h_operands);
    }
    //
    // BODY
    //
    make_round();
    // Filter out event queue for events matching the memory type.
    std::vector<event_t> sub_event_queue(event_queue);
    for (auto it = sub_event_queue.begin(); it != sub_event_queue.end(); ) {
        if (it->is_memory_x != is_memory_x) {
            it = sub_event_queue.erase(it);
        } else {
            it++;
        }
    }

    const size_t n_mt = mctr;
    const size_t n_et = sub_event_queue.size();
    
    int64_t event_ctr = 0,
            meas_ctr_offset = 0;
    for (size_t r = 0; r < rounds; r++) {
        push_back_range(program, round_program);
        // Make detection events.
        for (size_t i = 0; i < sub_event_queue.size(); i++) {
            event_t e = sub_event_queue.at(i);

            std::vector<int64_t> operands{event_ctr};
            for (int64_t m : e.m_operands) {
                int64_t mt = m + meas_ctr_offset;
                operands.push_back(mt);
                if (e.is_cross_round && r > 0) {
                    operands.push_back(mt - n_mt);
                }
            }
            if (operands.empty()) continue;
            qes::Instruction<> inst("event", operands);
            // Add annotations and properties
            if (e.owner != nullptr) {
                inst.put(print_v(e.owner));
            }
            if (e.owner != nullptr && e.owner->is_check()) {
                if (network->check_color_map.count(e.owner)) {
                    inst.put("color", network->check_color_map.at(e.owner));
                }
                int64_t base = rounds > 1 ? i+n_et : i;
                inst.put("base", base);
            } else {
                inst.put("flag");
            }
            program.push_back(inst);
            event_ctr++;
        }
        meas_ctr_offset += n_mt;
    }
    //
    // EPILOGUE
    //
    if (is_memory_x) {
        safe_emplace_back(program, "h", h_operands);
    }
    safe_emplace_back(program, "measure", m_operands);

    const size_t prev_meas_ctr_offset = meas_ctr_offset - n_mt;
    for (size_t i = 0; i < sub_event_queue.size(); i++) {
        event_t e = sub_event_queue.at(i);
        if (e.owner == nullptr || !e.owner->is_check()) continue;

        std::vector<int64_t> operands{event_ctr};
        for (int64_t m : e.m_operands) {
            int64_t mt = m + prev_meas_ctr_offset;
            operands.push_back(mt);
        }
        // Add data qubit measurements.
        sptr<tanner::vertex_t> tv = raw_network->v_tanner_raw_map.at(e.owner);
        for (sptr<tanner::vertex_t> tx : tanner_graph->get_neighbors(tv)) {
            sptr<raw_vertex_t> rx = raw_network->v_tanner_raw_map.at(tx);
            operands.push_back(data_meas_order_map.at(rx) + meas_ctr_offset);
        }
        qes::Instruction<> inst("event", operands);
        // Add annotations and properties.
        inst.put(print_v(e.owner));
        if (network->check_color_map.count(e.owner)) {
            inst.put("color", network->check_color_map.at(e.owner));
        }
        int64_t base = rounds > 1 ? i+n_et : i;
        inst.put("base", base);

        program.push_back(inst);
        event_ctr++;
    }
    auto obs_list = tanner_graph->get_obs(is_memory_x);
    int64_t obs_ctr = 0;
    for (auto& obs : obs_list) {
        std::vector<int64_t> operands{obs_ctr};
        for (sptr<tanner::vertex_t> tx : obs) {
            sptr<raw_vertex_t> rx = raw_network->v_tanner_raw_map.at(tx);
            operands.push_back(data_meas_order_map.at(rx) + meas_ctr_offset);
        }
        safe_emplace_back(program, "obs", operands);
        obs_ctr++;
    }
    return program;
}

void
Scheduler::make_round() {
    if (round_has_been_generated) return;
    qes::Program<> program;
    std::set<sptr<raw_vertex_t>> finished_checks;

    // Reset all tracking values.
    cycle = 0;
    event_queue.clear();
    mctr = 0;
    meas_ctr_map.clear();
    active_role_map.clear();

#ifdef PROTEAN_PERF
    Timer timer;
    fp_t t;
#endif
    while (true) {
#ifdef PROTEAN_PERF
        std::cout << "============= CYCLE " << cycle << " ===============" << std::endl;
#endif
        checks_this_cycle.clear();
        std::set<int64_t> r_operands;
        std::set<sptr<raw_vertex_t>> roles_this_cycle;
        for (sptr<raw_vertex_t> rv : raw_network->get_vertices()) {
            if (rv->is_check()) {
                sptr<phys_vertex_t> pv = network->role_to_phys.at(rv);
                size_t cyc = pv->cycle_role_map.at(rv);
                if (cyc <= cycle && !finished_checks.count(rv) && !roles_this_cycle.count(rv)) {
                    checks_this_cycle.push_back(rv);
                    insert_range(roles_this_cycle, raw_network->get_support(rv).all);
                    // Reset all parity, flags, and proxies associated with this qubit.
                    r_operands.insert(qu(rv));
                    auto& support = raw_network->get_support(rv);
                    for (sptr<raw_vertex_t> rx : support.flags) {
                        r_operands.insert(qu(rx));
                    }
                }
            }
        }
        safe_emplace_back(program, "reset", r_operands);
        if (checks_this_cycle.empty()) break;
#ifdef PROTEAN_PERF
        timer.clk_start();
#endif
        build_preparation(program);
#ifdef PROTEAN_PERF
        t = timer.clk_end();
        std::cout << "[ scheduler ] cycle " << cycle << " | preparation took " << t*1e-9 << "s" << std::endl;
        timer.clk_start();
#endif
        build_body(program);
#ifdef PROTEAN_PERF
        t = timer.clk_end();
        std::cout << "[ scheduler ] cycle " << cycle << " | body took " << t*1e-9 << "s" << std::endl;
        timer.clk_start();
#endif
        build_teardown(program);
#ifdef PROTEAN_PERF
        t = timer.clk_end();
        std::cout << "[ scheduler ] cycle " << cycle << " | teardown took " << t*1e-9 << "s" << std::endl;
        timer.clk_start();
#endif
        std::set<sptr<raw_vertex_t>> qubits_with_events;
        for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
            finished_checks.insert(rpq);
            qubits_with_events.insert(rpq);
            insert_range(qubits_with_events, raw_network->get_support(rpq).flags);
        }
        for (sptr<raw_vertex_t> rv : qubits_with_events) {
            declare_event_for_qubit(rv);
        }
        cycle++;
    }
    // Remove redundant instructions.
    remove_redundant_gates(program);
    
    program[0].put("timing_error");
    round_program = std::move(program);
    round_has_been_generated = true;
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
        for (const auto& [ rsq, sch ] : cx_sch_map) {
            // Make sure rsq is valid.
            test_and_set_exit_on_fail(rsq, "build_body");
            // Now get the k-th data qubit.
            if (k >= sch.data_order.size()) continue;
            sptr<raw_vertex_t> rdq = sch.data_order.at(k);
            if (rdq == nullptr) continue;
            cx_arr.push_back({rdq, rsq, sch.check->qubit_type == raw_vertex_t::type::xparity});
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
            if (std::find(m_operands.begin(), m_operands.end(), qu(rfq)) == m_operands.end()) {
                push_back_measurement(m_operands, rfq);
                r_operands.push_back(qu(rfq));
                release_qubit(rfq);
            }
        }
    }
    safe_emplace_back(program, "measure", m_operands);
}

void
Scheduler::remove_redundant_gates(qes::Program<>& program) {
    // We track two things:
    //  (1) last_inst_index_map: the last instruction for 
    //      each qubit operand (index into the program).
    //      This is used to identify redundant gates.
    //  (2) cx_pair_map: the corresponding cx operand for
    //      a qubit operand for some instruction. This is
    //      used to detect redundant CX gates.
    //  (3) deleted_operand_map: tracks which operands to
    //      delete from the corresponding instructions.
    std::map<int64_t, size_t> last_inst_index_map;
    TwoLevelMap<size_t, int64_t, int64_t> cx_pair_map;
    std::map<size_t, std::set<int64_t>> deleted_operand_map;

    size_t n_red = 0;
    for (size_t i = 0; i < program.size(); i++) {
        auto& inst = program[i];
        if (inst.get_name() == "cx") {
            for (size_t ii = 0; ii < inst.get_number_of_operands(); ii += 2) {
                int64_t q1 = inst.get<int64_t>(ii),
                        q2 = inst.get<int64_t>(ii+1);
                if (last_inst_index_map.count(q1) && last_inst_index_map.count(q2)
                    && last_inst_index_map.at(q1) == last_inst_index_map.at(q2)) 
                {
                    size_t j = last_inst_index_map[q1];
                    auto& _inst = program[j];
                    if (_inst.get_name() == "cx" && cx_pair_map[j].count(q1) && cx_pair_map[j][q1] == q2) {
                        // Then, this is redundant.
                        insert_all(deleted_operand_map[i], {q1, q2});
                        insert_all(deleted_operand_map[j], {q1, q2});
                        // Erase last_inst_index_map entry.
                        last_inst_index_map.erase(q1);
                        last_inst_index_map.erase(q2);
                        n_red++;
                        continue;
                    }
                }
                // Otherwise, update tracking structures.
                last_inst_index_map[q1] = i;
                last_inst_index_map[q2] = i;
                cx_pair_map[i][q1] = q2;
            }
        } else {
            for (size_t ii = 0; ii < inst.get_number_of_operands(); ii++) {
                int64_t q = inst.get<int64_t>(ii);
                if (last_inst_index_map.count(q)) {
                    size_t j = last_inst_index_map[q];
                    auto& _inst = program[j];
                    if (inst.get_name() == _inst.get_name()) {
                        deleted_operand_map[i].insert(q);
                        deleted_operand_map[j].insert(q);
                        last_inst_index_map.erase(q);
                        n_red++;
                        continue;
                    }
                }
                last_inst_index_map[q] = i;
            }
        }
    }
    // Update instructions now.
    for (auto& [i, deleted] : deleted_operand_map) {
        auto& inst = program[i];
        std::vector<int64_t> new_operands;
        for (size_t ii = 0; ii < inst.get_number_of_operands(); ii++) {
            int64_t q = inst.get<int64_t>(ii);
            if (!deleted.count(q)) new_operands.push_back(q);
        }
        qes::Instruction<> _inst(inst.get_name(), new_operands);
        // Copy over annotations and properties.
        for (std::string a : inst.get_annotations()) _inst.put(a);
        for (auto& [k, v] : inst.get_property_map()) _inst.put(k, v);
        inst = std::move(_inst);
    }
    // Finally, do a final pass and remove any bad instructions.
    for (auto it = program.begin(); it != program.end(); ) {
        if (it->get_number_of_operands() == 0) it = program.erase(it);
        else it++;
    }
    std::cout << "[ status ] deleted " << n_red << " redundant instructions" << std::endl;
}

void
Scheduler::prep_tear_h_gates(qes::Program<>& program) {
    std::set<int64_t> h_operands;
    for (sptr<raw_vertex_t> rpq : checks_this_cycle) {
        test_and_set_exit_on_fail(rpq, "build_preparation");
        // Perform an H gate if necessary.
        auto& support = get_support(rpq);
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            h_operands.insert(qu(rpq));
        } else {
            for (sptr<raw_vertex_t> rfq : support.flags) {
                test_and_set_exit_on_fail(rfq, "build_preparation");
                h_operands.insert(qu(rfq));
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
            cx_arr.push_back({rfq, rpq, rpq->qubit_type == raw_vertex_t::type::xparity});
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
        sptr<raw_vertex_t> src = cx.src,
                                dst = cx.dst;
        path_arr[i] = raw_network->get_proxy_walk_path(src, dst);
    }
    // Also track undo paths:
    std::vector<std::vector<sptr<raw_vertex_t>>> undo_arr;
    std::vector<bool> undo_is_for_x_check;
    std::vector<size_t> undo_k_arr;
    // Now perform CX along these paths.
    std::set<sptr<raw_vertex_t>> proxies_in_use;
    while (true) {
        std::set<int64_t> in_use;
        std::vector<int64_t> h_pre_operands, h_post_operands, cx_operands, r_operands;
        for (size_t i = 0; i < cx_arr.size(); i++) {
            bool is_for_x_check = cx_arr.at(i).is_for_x_check;
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
                    },
                    is_for_x_check)
                )
            {
                continue;
            }

            sptr<raw_vertex_t> ry = path.at(k);
            if (k < path.size()-1 && ry->qubit_type == raw_vertex_t::type::proxy) {
                proxies_in_use.insert(ry);
                // If this is for an x check, then perform an H gate on path[k] as well.
                if (is_for_x_check) {
                    h_pre_operands.push_back(qu(ry));
                }
            }
            // Check if we have a valid undo path. If so, mark it as a path
            // we need to complete.
            k++;
            if (k == path.size() && path.size() > 2) {
                undo_arr.emplace_back(path.rbegin()+1, path.rend());
                undo_k_arr.push_back(1);
                undo_is_for_x_check.push_back(is_for_x_check);
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
                    !undo_is_for_x_check[i])
                )
            {
                continue;
            }

            sptr<raw_vertex_t> rx = path.at(k-1);
            if (rx->qubit_type == raw_vertex_t::type::proxy) {
                proxies_in_use.erase(rx);
                if (undo_is_for_x_check[i]) {
                    h_post_operands.push_back(qu(rx));
                }
            }
            k++;
        }
        if (cx_operands.empty()) break;
        safe_emplace_back(program, "h", h_pre_operands);
        safe_emplace_back(program, "cx", cx_operands);
        safe_emplace_back(program, "h", h_post_operands);
    }
}

std::map<sptr<raw_vertex_t>, sch_data_t>
Scheduler::compute_schedules() {
    std::map<sptr<raw_vertex_t>, sch_data_t> sch_map;
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
    TwoLevelMap<sptr<raw_vertex_t>, sptr<raw_vertex_t>, size_t>
        sch_time_map;

    auto& fam = raw_network->flag_assignment_map;
    std::vector<sptr<raw_vertex_t>> sorted_checks(checks_this_cycle);
    std::sort(sorted_checks.begin(), sorted_checks.end(),
            [&] (sptr<raw_vertex_t> rx, sptr<raw_vertex_t> ry)
            {
                return raw_network->get_support(rx).data.size()
                        < raw_network->get_support(ry).data.size();
            });
#ifdef PROTEAN_PERF
    Timer timer;
    timer.clk_start();
#endif
    for (sptr<raw_vertex_t> rpq : sorted_checks) {
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
            sptr<raw_vertex_t> rx_rsq = fam[rpq].count(rx) ? fam[rpq][rx] : rpq;
            if (sch_map.count(rx_rsq)) {
                // Then rx has been scheduled already. We can skip adding constraints
                // for this variable.
                lp_constr_t con(x, sch_time_map[rx_rsq][rx], lp_constr_t::direction::eq);
                LP.add_constraint(con);
                continue;
            }
            for (sptr<raw_vertex_t> ry : support.data) {
                if (rx >= ry) continue;
                lp_var_t y = LP.get_var(ry);
                sptr<raw_vertex_t> ry_rsq = fam[rpq].count(ry) ? fam[rpq][ry] : rpq;
                // Only add the following constraint if rx and ry both interact
                // with rpq or the same flag.
                if (rx_rsq == ry_rsq) {
                    lp_constr_t con(x, y, lp_constr_t::direction::neq);
                    LP.add_constraint(con);
                }
            }
        }
        // Add stabilizer constraints:
        // We will require commutation on a flag granularity.
        TwoLevelMap<sptr<raw_vertex_t>, sptr<raw_vertex_t>, lp_expr_t> ind_sum_map;
        TwoLevelMap<sptr<raw_vertex_t>, sptr<raw_vertex_t>, size_t> ind_ctr_map;
        for (sptr<raw_vertex_t> rdq : support.data) {
            lp_var_t x = LP.get_var(rdq);

            sptr<raw_vertex_t> rsq = fam[rpq].count(rdq) ? fam[rpq][rdq] : rpq;
            for (auto& tt : existing_cnot_times[rdq]) {
                size_t t = std::get<0>(tt);
                sptr<raw_vertex_t> owner = std::get<1>(tt);
                sptr<raw_vertex_t> _rsq = fam[owner].count(rdq) ? fam[owner][rdq] : owner;
                // Add conflict constraint if owner and rpq are not using the same flag.
                if (rsq != _rsq) {
                    lp_constr_t con1(x, t, lp_constr_t::direction::neq);
                    LP.add_constraint(con1);
                }
                // if is_x_check != tt_from_x, then we need to add commutation
                // constraints.
                if (rpq->qubit_type == owner->qubit_type) continue;
                lp_expr_t& ind_sum = ind_sum_map[_rsq][rsq];

                const double M = 100'000;
                lp_var_t y = LP.add_slack_var(0, 1, lp_var_t::bounds::both, lp_var_t::domain::binary);
                // x - t <= My and t - x <= M(1-y).
                lp_constr_t con2(x - t, M*y, lp_constr_t::direction::le),
                            con3(t - x, M*(1-y), lp_constr_t::direction::le);
                LP.add_constraint(con2);
                LP.add_constraint(con3);

                ind_sum += y;
                ind_ctr_map[_rsq][rsq]++;
                // We have added the constraint for _rsq and rsq. But now we need to add it
                // for owner and rpq if _rsq != owner or rsq != rpq.
                if (_rsq != owner) {
                    ind_sum_map[owner][rsq] += y;
                    ind_ctr_map[owner][rsq]++;
                    if (rsq != rpq) {
                        ind_sum_map[owner][rpq] += y;
                        ind_ctr_map[owner][rpq]++;
                    }
                }
                if (rsq != rpq) {
                    ind_sum_map[_rsq][rpq] += y;
                    ind_ctr_map[_rsq][rpq]++;
                }
            }
        }
        for (auto& p1 : ind_sum_map) {
            for (auto& p2 : p1.second) {
                lp_expr_t& ind_sum = p2.second;
                if (ind_ctr_map[p1.first][p2.first] & 1) {
                    continue;
                }
                // Constrain ind_sum to be even.
                lp_var_t y = LP.add_slack_var(0, 0, lp_var_t::bounds::lower, lp_var_t::domain::integer);
                lp_constr_t con(ind_sum, 2*y, lp_constr_t::direction::eq);
                LP.add_constraint(con);
            }
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
                if (fam[rpq].count(rdq)) {
                    sptr<raw_vertex_t> rsq = fam[rpq][rdq];
                    if (sch_map.count(rsq)) {
                        std::cerr << " | fixed to be " << sch_time_map[rsq][rdq];
                    }
                    std::cerr << " | rsq = " << print_v(rsq);
                }
                std::cerr << std::endl;
            }
            std::cerr << "possible problems:" << std::endl
                    << "\tnumber of variables (" << LP.variables.size() << ")" << std::endl
                    << "\tnumber of constraints (" << LP.l_constraints.size() << ")" << std::endl;
            exit(1);
        }

        for (sptr<raw_vertex_t> rdq : support.data) {
            sptr<raw_vertex_t> rsq = fam[rpq].count(rdq) ? fam[rpq][rdq] : rpq;
            size_t t = static_cast<size_t>(round(LP.get_value(rdq)));
            existing_cnot_times[rdq].emplace_back(t, rpq);

            if (!sch_map.count(rsq)) {
                sch_map[rsq] = { rpq, std::vector<sptr<raw_vertex_t>>(upper_bound, nullptr) };
            }
            sch_map[rsq].data_order[t-1] = rdq;
            tlm_put(sch_time_map, rsq, rdq, t);
        }
    }
#ifdef PROTEAN_PERF
    fp_t t = timer.clk_end();
    std::cout << "[ compute_schedules ] computing schedule took " << t*1e-9 << "s" << std::endl;
#endif
    return sch_map;
}

}   // protean
}   // qontra
