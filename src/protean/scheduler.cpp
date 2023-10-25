/*
 *  author: Suhas Vittal
 *  date:   23 October 2023
 * */

#include "protean/scheduler.h"

namespace qontra {

using namespace graph;

namespace protean {

css_code_data_t
compute_schedule_from_tanner_graph(TannerGraph& tanner_graph) {
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
                    e->src = (void*) vj;
                    e->is_undirected = true;
                    gr.add_edge(e);
                }
            }
        }
    }
    // Now compute the schedule via a BFS on the graph.
    std::deque<stab_vertex_t*> bfs{vertices[0]};
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
            int t = (int) mgr->get_value(q);
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
        lp_constr_t<uint> con(max_of_all, q, ConstraintDirection::ge);
        mgr->add_constraint(con);
    }
    // Add uniqueness constraint (any qi != qj).
    const uint w = support.size();
    for (uint i = 0; i < w; i++) {
        lp_var_t<uint>* q1v = mgr->get_var(support[i]);
        for (uint j = i+1; j < w; j++) {
            lp_var_t<uint>* q2v = mgr->get_var(support[j]);
            lp_constr_t<uint> con(q1v, q2v, ConstraintDirection::neq);
            mgr->add_constraint(con);
        }
    }

    // Setup constraints with other stabilizers.
    for (auto s : gr.get_neighbors(s0)) {
        const auto& other_support = s->support;
        const auto& sch = s->sch_qubit_to_time;
        // We only care about neighboring stabilizers if they
        // already have a schedule.
        if (sch.empty()) continue;
        bool sum_is_nonzero = false;
        lp_expr_t<uint> ind_sum;
        for (auto q : support) {
            if (std::find(other_support.begin(), other_support.end(), q) == other_support.end()) {
                continue;
            }
            lp_var_t<uint>* qv = mgr->get_var(q);
            int t = sch.at(q);
            lp_constr_t<uint> con1(qv, t, ConstraintDirection::neq);
            mgr->add_constraint(con1);
            if (s0->qubit_to_pauli[q] != s->qubit_to_pauli[q]) {
                const double M = 100000;
                lp_var_t<uint>* ind = mgr->add_slack_var(0, 1, VarBounds::both, VarDomain::binary);
                lp_constr_t<uint> con2(qv - t, -M*lp_expr_t<uint>(ind), ConstraintDirection::le);
                mgr->add_constraint(con2);
                // Also add ind to ind_sum.
                ind_sum += ind;
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
