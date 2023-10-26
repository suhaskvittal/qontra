/* author: Suhas Vittal date:   23 October 2023
 * */

#include "protean/scheduler.h"

namespace qontra {

using namespace graph;
using namespace enumerator;

namespace protean {

mq_pauli_op_t
mul(const mq_pauli_op_t& p1, const mq_pauli_op_t& p2) {
    std::set<pauli_op_t> ss(p1.begin(), p1.end());
    for (pauli_op_t x : p2) xor_into(ss, x);
    return mq_pauli_op_t(ss.begin(), ss.end());
}

int
count_anticommutations(const mq_pauli_op_t& p1, const mq_pauli_op_t& p2) {
    std::map<uint, pauli> qubit_to_pauli;
    int c = 0;
    for (pauli_op_t x : p1) qubit_to_pauli[x.second] = x.first;
    for (pauli_op_t x : p2) {
        if (qubit_to_pauli.count(x.second) && qubit_to_pauli[x.second] != x.first) {
            c++;
        }
    }
    return c;
}

std::vector<stabilizer_t>
get_stabilizers(css_code_data_t code_data) {
    std::vector<stabilizer_t> stabs;
    for (uint pq : code_data.xparity_list) {
        stabilizer_t s;
        for (int32_t dq : code_data.check_schedules[pq]) {
            if (dq < 0) continue;
            s.push_back(std::make_pair(pauli::x, (uint)dq));
        }
        stabs.push_back(s);
    }

    for (uint pq : code_data.zparity_list) {
        stabilizer_t s;
        for (int32_t dq : code_data.check_schedules[pq]) {
            if (dq < 0) continue;
            s.push_back(std::make_pair(pauli::z, (uint)dq));
        }
        stabs.push_back(s);
    }
    return stabs;
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
    // Set the observables
    for (auto obs : tanner_graph.get_obs(true)) {
        std::vector<uint> qlist;
        for (auto dv : obs) qlist.push_back(tanner_vertex_to_qubit[dv]);
        code_data.x_obs_list.push_back(qlist);
    }
    for (auto obs : tanner_graph.get_obs(false)) {
        std::vector<uint> qlist;
        for (auto dv : obs) qlist.push_back(tanner_vertex_to_qubit[dv]);
        code_data.z_obs_list.push_back(qlist);
    }
    // Flag qubits will remain empty for now.
    return code_data;
}

css_code_data_t
make_fault_tolerant(css_code_data_t code_data) {
    bool _g_use_mpi = experiments::G_USE_MPI;
    experiments::G_USE_MPI = false;

    std::vector<stabilizer_t> stab = get_stabilizers(code_data);
    // Finally, when introducing flag qubits for a check C and data qubit
    // D1, we will connect the flag to the upcoming data qubit D2 in the
    // schedule. Hence, we maintain a next qubit table.
    TwoLevelMap<uint, uint, uint> next_qubit_in_schedule;
    for (auto& check_sch_pair : code_data.check_schedules) {
        uint check = check_sch_pair.first;
        auto& data_qubit_order = check_sch_pair.second;
        int32_t prev = -1;
        for (int32_t dq : data_qubit_order) {
            if (dq >= 0) {
                if (prev >= 0) {
                    tlm::put(next_qubit_in_schedule, check, (uint)prev, (uint)dq);
                }
                prev = dq;
            }
        }
    }
    // We will analyze the errors on the syndrome extraction schedule.
    schedule_t sxt_sch = write_syndrome_extraction_ops(code_data);
    std::vector<error_record_t> all_errors = enumerate_errors(sxt_sch);
    // Now, we will examine these errors and introduce flag qubits wherever
    // necessary.
    uint flag_qubit_ctr = code_data.data_qubits.size() + code_data.parity_qubits.size();
    for (error_record_t& rc : all_errors) {
        Instruction& inst = std::get<1>(rc);
        uint qubit = (uint)std::get<2>(rc);
        if (inst.name != "cx") continue; // Nothing we can do.
        stim::simd_bits x_errors = std::get<4>(rc);
        stim::simd_bits z_errors = std::get<5>(rc);
        // Compute pauli operator for error.
        mq_pauli_op_t error_op;
        for (uint dq : code_data.data_qubits) {
            if (x_errors[dq] & z_errors[dq])    error_op.push_back(std::make_pair(pauli::y, dq));
            else if (x_errors[dq])              error_op.push_back(std::make_pair(pauli::x, dq));
            else if (z_errors[dq])              error_op.push_back(std::make_pair(pauli::z, dq));
        }
        if (error_op.size() == 1) continue; // Not an issue.
        // Check if the error operator forms a stabilizer generator.
        // Ideally, we'd like to check if it forms a stabilizer, but that's
        // a bit harder.
        //
        // Also compute minimum weight dual error (product of stab with error).
        //
        // We can exit early if w == 1 (not an issue) or w == 0 (error is a
        // stabilizer).
        mq_pauli_op_t min_weight_dual;
        uint w = std::numeric_limits<uint>::max();
        for (const stabilizer_t& s : stab) {
            mq_pauli_op_t dual = mul(error_op, s);
            if (dual.size() < w) {
                min_weight_dual = dual;
                w = dual.size();
                if (w <= 1) break;
            }
        }
        if (w <= 1) continue;
        if (min_weight_dual.size() < error_op.size()) error_op = min_weight_dual;
        // Now, we must analyze and check for errors. Check against the
        // observables.
        pauli pauli_array[] = { pauli::x, pauli::z };
        std::vector<std::vector<uint>> obs_list_array[] = { code_data.x_obs_list, code_data.z_obs_list };
        bool is_fault_tolerant = true;
        for (int i = 0; i < 2; i++) {
            pauli p = pauli_array[i];
            const auto& obs_list = obs_list_array[i];
            for (const auto& obs : obs_list) {
                mq_pauli_op_t obs_op;
                for (uint q : obs) obs_op.push_back(std::make_pair(p, q));
                int c = count_anticommutations(error_op, obs_op);
                if (c > 1) {
                    is_fault_tolerant = false;
                    goto obs_cmp_exit;
                }
            }
        }
obs_cmp_exit:
        if (is_fault_tolerant) continue;
        // Check CX instruction
        uint other_qubit;
        // Find other other qubit in the CNOT.
        for (uint i = 0; i < inst.operands.qubits.size(); i += 2) {
            uint q1 = inst.operands.qubits[i];
            uint q2 = inst.operands.qubits[i+1];
            if (q1 == qubit) {
                other_qubit = q2;
                break;
            }
            if (q2 == qubit) {
                other_qubit = q1;
                break;
            }
        }
        // One of the qubit operands is a check qubit.
        uint data_qubit = qubit, check_qubit = other_qubit;
        if (code_data.check_schedules.count(qubit)) {
            data_qubit = other_qubit;
            check_qubit = qubit;
        }
        // If the data qubit is already assigned a flag qubit, move on.
        if (code_data.flag_usage[check_qubit].count(data_qubit)) continue;
        uint next_data_qubit = next_qubit_in_schedule[check_qubit][data_qubit];
        // It should never be the case that the next data qubit is already
        // assigned a flag.
        tlm::put(code_data.flag_usage, check_qubit, data_qubit, flag_qubit_ctr);
        tlm::put(code_data.flag_usage, check_qubit, next_data_qubit, flag_qubit_ctr);
        if (std::find(code_data.xparity_list.begin(), code_data.xparity_list.end(), check_qubit)
                == code_data.xparity_list.end())
        {
            code_data.xflag_list.push_back(flag_qubit_ctr);
        } else {
            code_data.zflag_list.push_back(flag_qubit_ctr);
        }
        code_data.parity_to_flags[check_qubit].push_back(flag_qubit_ctr);
        code_data.flag_qubits.push_back(flag_qubit_ctr);
        flag_qubit_ctr++;
    }

    experiments::G_USE_MPI = _g_use_mpi;
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
