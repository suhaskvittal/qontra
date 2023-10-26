/*
 *  author: Suhas Vittal
 *  date:   23 October 2023
 * */

#include "protean/scheduler.h"

namespace qontra {

using namespace graph;

namespace protean {

schedule_t
write_memory_experiment(css_code_data_t code_data, uint rounds, bool is_memory_x) {
    schedule_t prog;
    //
    // PROLOGUE.
    //
    std::vector<uint> all_qubits(code_data.data_qubits);
    for (auto q : code_data.parity_qubits) all_qubits.push_back(q);
    for (auto q : code_data.flag_qubits) all_qubits.push_back(q);

    prog.push_back(Instruction::gate("reset", all_qubits));
    if (is_memory_x) {
        prog.push_back(Instruction::gate("h", code_data.data_qubits));
    }
    //
    // BODY
    //
    std::vector<uint> hlist(code_data.xparity_list);
    for (auto q : code_data.zflag_list) hlist.push_back(q);

    // These are the H gates used at the start and end of the round.
    Instruction hstart = Instruction::gate("h", hlist);
    hstart.annotations.insert("inject_timing_error");
    Instruction hend = Instruction::gate("h", hlist);

    // This is the set of CX gates required to entangle/disentangle the flags
    // from the parity qubits.
    std::vector<Instruction> flagcx_array;
    int index = 0;
    while (true) {
        std::vector<uint> qarr;
        bool had_flag = false;
        for (auto pq : code_data.zparity_list) {
            std::vector<uint> flags = code_data.parity_to_flags[pq];
            if (index >= flags.size())  continue;
            had_flag = true;

            uint fq = flags[index];
            qarr.push_back(fq);
            qarr.push_back(pq);
        }
        for (auto pq : code_data.xparity_list) {
            std::vector<uint> flags = code_data.parity_to_flags[pq];
            if (index >= flags.size())  continue;
            had_flag = true;

            uint fq = flags[index];
            qarr.push_back(pq);
            qarr.push_back(fq);
        }
        if (!had_flag) break;
        flagcx_array.push_back(Instruction::gate("cx", qarr));
        index++;
    }

    // Finally, the proper syndrome extraction CNOTs.
    std::vector<Instruction> extcx_array;
    index = 0;
    while (true) {
        std::vector<uint> qarr;
        bool had_op = false;
        std::cout << "index = " << index << "\n";
        for (auto pq : code_data.zparity_list) {
            std::vector<int32_t> supp = code_data.check_schedules[pq];
            if (index >= supp.size())   continue;
            had_op = true;
            if (supp[index] < 0) continue;

            uint dq = supp[index];
            // This may have a flag qubit -- check if it does.
            if (code_data.flag_usage[pq].count(dq)) {
                uint fq = code_data.flag_usage[pq][dq];
                qarr.push_back(dq);
                qarr.push_back(fq);
            } else {
                qarr.push_back(dq);
                qarr.push_back(pq);
            }
        }

        for (auto pq : code_data.xparity_list) {
            std::vector<int32_t> supp = code_data.check_schedules[pq];
            if (index >= supp.size())   continue;
            had_op = true;
            if (supp[index] < 0) continue;

            uint dq = supp[index];
            // This may have a flag qubit -- check if it does.
            if (code_data.flag_usage[pq].count(dq)) {
                uint fq = code_data.flag_usage[pq][dq];
                qarr.push_back(fq);
                qarr.push_back(dq);
            } else {
                qarr.push_back(pq);
                qarr.push_back(dq);
            }
        }
        if (!had_op) break;
        extcx_array.push_back(Instruction::gate("cx", qarr));
        index++;
    }
    std::vector<uint> all_flag_parity_qubits;
    std::map<uint, uint> qubit_to_meas_time;
    uint64_t meas_per_round = 0;
    for (uint q : code_data.zparity_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.xparity_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.zflag_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.xflag_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    Instruction flagparitymeas = Instruction::gate("measure", all_flag_parity_qubits);
    Instruction flagparityreset = Instruction::gate("reset", all_flag_parity_qubits);
    // Now add all the instructions to the schedule.
    std::vector<uint> event_qubits;
    if (is_memory_x) {
        for (uint q : code_data.xparity_list) event_qubits.push_back(q);
        for (uint q : code_data.zflag_list) event_qubits.push_back(q);
    } else {
        for (uint q : code_data.zparity_list) event_qubits.push_back(q);
        for (uint q : code_data.xflag_list) event_qubits.push_back(q);
    }
    uint64_t event_ctr = 0;
    for (uint r = 0; r < rounds; r++) {
        prog.push_back(hstart);
        for (auto cx_inst : flagcx_array)   prog.push_back(cx_inst);
        for (auto cx_inst : extcx_array)    prog.push_back(cx_inst);
        for (auto cx_inst : flagcx_array)   prog.push_back(cx_inst);
        prog.push_back(hend);
        prog.push_back(flagparitymeas);
        prog.push_back(flagparityreset);
        // Create detection events.
        for (uint i = 0; i < event_qubits.size(); i++) {
            uint q = event_qubits[i];
            uint mt = qubit_to_meas_time[q];
            if (r > 0) {
                prog.push_back(Instruction::event(event_ctr++, {mt + r*meas_per_round, mt + (r-1)*meas_per_round}));
            } else {
                prog.push_back(Instruction::event(event_ctr++, {mt}));
            }
        }
    }
    const uint64_t total_measurements = rounds*meas_per_round;
    //
    // Epilogue
    //
    if (is_memory_x) {
        prog.push_back(Instruction::gate("h", code_data.data_qubits));
    }
    // We need to record when we measure the data qubits when creating the final
    // detection events and the observables.
    for (uint i = 0; i < code_data.data_qubits.size(); i++) {
        uint q = code_data.data_qubits[i];
        qubit_to_meas_time[q] = total_measurements+i;
    }
    prog.push_back(Instruction::gate("measure", code_data.data_qubits));
    // Create final detection events.
    std::vector<uint> relevant_checks;
    if (is_memory_x) relevant_checks = code_data.xparity_list;
    else             relevant_checks = code_data.zparity_list;
    for (uint pq : relevant_checks) {
        std::vector<int32_t> supp = code_data.check_schedules[pq];
        std::vector<uint> mtimes;
        for (int32_t dq : supp) {
            if (dq < 0) continue;
            uint mt = qubit_to_meas_time[dq];
        }
        mtimes.push_back((rounds-1)*meas_per_round + qubit_to_meas_time[pq]);
        prog.push_back(Instruction::event(event_ctr++, mtimes));
    }
    // Create logical observables.
    std::vector<std::vector<uint>> obs_list;
    if (is_memory_x)    obs_list = code_data.x_obs_list;
    else                obs_list = code_data.z_obs_list;
    uint64_t obs_ctr = 0;
    for (auto obs : obs_list) {
        std::vector<uint> mtimes;
        for (uint dq : obs) mtimes.push_back(qubit_to_meas_time[dq]);
        prog.push_back(Instruction::obs(obs_ctr++, mtimes));
    }
    // The schedule is done -- we can return now.
    return prog;
}

css_code_data_t
compute_schedule_from_tanner_graph(TannerGraph& tanner_graph, int start) {
    // Translate tanner graph to a StabilizerGraph
    StabilizerGraph gr;
    int max_stab_weight = 0;
    for (auto pv : tanner_graph.get_checks()) {
        stabilizer_t stab;
        support_t support;
        std::map<uint, pauli> q2p;
        for (auto dv : tanner_graph.get_neighbors(pv)) {
            pauli p = pauli::x;
            if (pv->qubit_type == tanner::vertex_t::Type::zparity) {
                p = pauli::z;
            }
            uint q = (uint) (tanner::VERTEX_ID_NUMBER_MASK & dv->id);
            pauli_op_t op = std::make_pair(p, q);
            stab.push_back(op);
            support.push_back(q);
            q2p[q] = p;
        }
        stab_vertex_t* v = new stab_vertex_t;
        v->id = pv->id;
        v->check = pv;
        v->stabilizer = stab;
        v->support = support;
        v->qubit_to_pauli = q2p;
        gr.add_vertex(v);

        if (support.size() > max_stab_weight) {
            max_stab_weight = support.size();
        }
    }
    // Add edges to the StabilizerGraph.
    auto vertices = gr.get_vertices();
    for (uint i = 0; i < vertices.size(); i++) {
        auto vi = vertices[i];
        for (uint j = i+1; j < vertices.size(); j++) {
            auto vj = vertices[j];
            // Check for intersections in the support.
            for (auto qi : vi->support) {
                pauli pi = vi->qubit_to_pauli[qi];
                if (!vj->qubit_to_pauli.count(qi))  continue;
                
                pauli pj = vj->qubit_to_pauli[qi];
                if (pi != pj) {
                    // Make edge -- the paulis anticommute.
                    stab_edge_t* e = new stab_edge_t;
                    e->src = (void*) vi;
                    e->dst = (void*) vj;
                    e->is_undirected = true;
                    gr.add_edge(e);
                    break;
                }
            }
        }
    }
    // Now compute the schedule via a BFS on the graph.
    std::deque<stab_vertex_t*> bfs{vertices[start]};
    std::set<stab_vertex_t*> visited;
    int max_time = 0;
    while (bfs.size()) {
        auto s = bfs.front();
        bfs.pop_front();
        if (visited.count(s))   continue;
        // Compute schedule for this stabilizer.
        // First construct an ILP.
        LPManager<uint>* mgr = construct_scheduling_program(s, gr, max_stab_weight);
        mgr->solve();
        // Get variable values.
        for (uint q : s->support) {
            int t = (int) round(mgr->get_value(q));
            s->sch_qubit_to_time[q] = t;
            s->sch_time_to_qubit[t] = q;
            if (t > max_time) max_time = t;
        }
        // Traverse to neighbors.
        for (auto u : gr.get_neighbors(s)) {
            bfs.push_back(u);
        }
        visited.insert(s);
        delete mgr;
    }
    // Now that we have computed all the schedules, let us compile the code
    // data.
    css_code_data_t code_data;
    code_data.schedule_depth = max_time;
    // We will maintain an assignment of tanner graph vertices to qubits.
    std::map<tanner::vertex_t*, uint64_t> tanner_vertex_to_qubit;
    uint64_t next_qubit_id;
    for (auto dv : tanner_graph.get_vertices_by_type(tanner::vertex_t::Type::data)) {
        code_data.data_qubits.push_back(next_qubit_id);
        tanner_vertex_to_qubit[dv] = next_qubit_id;
        next_qubit_id++;
    }
    for (auto pv : tanner_graph.get_checks()) {
        code_data.parity_qubits.push_back(next_qubit_id);
        tanner_vertex_to_qubit[pv] = next_qubit_id;
        if (pv->qubit_type == tanner::vertex_t::Type::xparity) {
            code_data.xparity_list.push_back(next_qubit_id);
        } else {
            code_data.zparity_list.push_back(next_qubit_id);
        }
        stab_vertex_t* sv = gr.get_vertex(pv->id);
        const auto& sch = sv->sch_time_to_qubit;
        for (int t = 1; t <= max_time; t++) {
            if (!sch.count(t)) {
                code_data.check_schedules[next_qubit_id].push_back(-1);
            } else {
                code_data.check_schedules[next_qubit_id].push_back(sch.at(t));
            }
        }
        next_qubit_id++;
    }
    // Flag qubits will remain empty for now.
    return code_data;
}

LPManager<uint>*
construct_scheduling_program(
        stab_vertex_t* s0,
        StabilizerGraph& gr,
        int max_stabilizer_weight)
{
    LPManager<uint>* mgr = new LPManager<uint>;

    // This variable is greater than all qubit variables.
    lp_var_t<uint>* max_of_all = mgr->add_slack_var(1, 2*max_stabilizer_weight, VarBounds::both, VarDomain::integer);

    const auto& support = s0->support;
    // Create variables for each qubit.
    for (uint q : support) {
        lp_var_t<uint>* qv = mgr->add_var(q, 1, 2*max_stabilizer_weight, VarBounds::both, VarDomain::integer);
        lp_constr_t<uint> con(lp_expr_t<uint>(max_of_all) - lp_expr_t<uint>(qv), 0.0, ConstraintDirection::ge);
        mgr->add_constraint(con);
    }
    // Add uniqueness constraint (any qi != qj).
    const uint w = support.size();
    for (uint i = 0; i < w; i++) {
        lp_var_t<uint>* q1v = mgr->get_var(support[i]);
        for (uint j = i+1; j < w; j++) {
            lp_var_t<uint>* q2v = mgr->get_var(support[j]);
            lp_constr_t<uint> con(lp_expr_t<uint>(q1v), lp_expr_t<uint>(q2v), ConstraintDirection::neq);
            mgr->add_constraint(con);
        }
    }

    // Setup constraints with other stabilizers.
    for (auto s : gr.get_vertices()) {
        if (s == s0) continue;
        const auto& other_support = s->support;
        const auto& sch = s->sch_qubit_to_time;
        // We only care about neighboring stabilizers if they
        // already have a schedule.
        if (sch.empty()) continue;
        bool sum_is_nonzero = false;
        lp_expr_t<uint> ind_sum;
        for (auto q : support) {
            if (!sch.count(q)) {
                continue;
            }
            lp_var_t<uint>* qv = mgr->get_var(q);
            int t = sch.at(q);
            lp_constr_t<uint> con1(qv, t, ConstraintDirection::neq);
            mgr->add_constraint(con1);
            if (s0->qubit_to_pauli[q] != s->qubit_to_pauli[q]) {
                const double M = 100000;
                lp_var_t<uint>* ind = mgr->add_slack_var(0, 1, VarBounds::both, VarDomain::binary);
                // x - t <= My
                lp_constr_t<uint> con2(
                        lp_expr_t<uint>(qv) - t, M*lp_expr_t<uint>(ind), ConstraintDirection::le);
                // t - x <= M(1-y)
                lp_constr_t<uint> con3(
                        t - lp_expr_t<uint>(qv), M - M*lp_expr_t<uint>(ind), ConstraintDirection::le);
                mgr->add_constraint(con2);
                mgr->add_constraint(con3);
                // Also add ind to ind_sum.
                ind_sum += lp_expr_t<uint>(ind);
                sum_is_nonzero = true;
            }
        }
        if (sum_is_nonzero) {
            lp_var_t<uint>* sum_ind_div_two = mgr->add_slack_var(0, 0, VarBounds::lower, VarDomain::integer);
            lp_constr_t<uint> con(ind_sum, 2*lp_expr_t<uint>(sum_ind_div_two), ConstraintDirection::eq);
            mgr->add_constraint(con);
        }
    }
    // Build the LP.
    mgr->build(max_of_all, false);
    return mgr;
}

}   // protean
}   // qontra
