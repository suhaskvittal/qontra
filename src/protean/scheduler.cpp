/*
 *  author: Suhas Vittal
 *  date:   23 October 2023
 * */

#ifndef PROTEAN_SCHEDULER_h
#define PROTEAN_SCHEDULER_h

#include "protean/scheduler.h"

namespace qontra {

using namespace graph;

namespace protean {

static IloEnv env;

void finalize_cplex() { env.end(); }

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
            if (pv->qubit_type == tanner::vertex_t::Type::zparity) p = pauli::z;
            uint q = (uint) (tanner::VERTEX_ID_NUMBER_MASK & dv->id);
            pauli_op_t operator = std::make_pair(p, q);
            stab.push_back(operator);
            support.push_back(q);
            p2q[q] = p;
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
    auto& vertices = gr.get_vertices();
    for (uint i = 0; i < vertices.size(); i++) {
        auto vi = vertices[i];
        for (uint j = i+1; j < vertices.size(); j++) {
            auto vj = vertices[j];
            // Check for intersections in the support.
            for (auto qi : vi->support) {
                pauli pi = vi->qubit_to_pauli[qi];
                if (vj->pauli_to_qubit.count(qi)) {
                    pauli pj = vj->pauli_to_qubit[qi];
                    if (pi != pj) {
                        // Make edge -- the paulis anticommute.
                        stab_edge_t e = new stab_edge_t;
                        e->src = (void*) vi;
                        e->src = (void*) vj;
                        e->is_undirected = true;
                        gr.add_edge(e);
                    }
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
        std::map<uint, IloIntVar> qubit_to_variable;
        IloModel prog = construct_scheduling_program(s, gr, max_stab_weight, qubit_to_variable);
        IloCplex solver(prog);
        if (!prog.solve()) {
            return (css_code_data_t) {}; // Infeasible solution found.
        }
        // Get variable values.
        for (uint q : s->support) {
            int t = (int) solver.getValue(qubit_to_variable[q]);
            s->sch_qubit_to_time[q] = t;
            s->sch_time_to_qubit[t] = q;
            if (t > max_time) max_time = t;
        }
        // Traverse to neighbors.
        for (auto u : gr.get_neighbors(s)) {
            bfs.push_back(u);
        }
        visited.insert(s);
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
        const auto& sch = sv->time_to_qubit;
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

IloModel
construct_scheduling_program(
        stab_vertex_t* s0,
        StabilizerGraph& gr,
        int max_stabilizer_weight,
        std::map<uint, IloIntVar>& qubit_to_variable)
{
    IloModel model(env);

    const auto& support = s0->support;
    IloIntVar max_of_all(env, 1, 2*max_stabilizer_weight);
    model.add(max_of_all);
    model.add(IloMinimize(max_of_all));

    IloIntVarArray all_qubit_vars(env);

    for (uint q : support) {
        IloIntVar x(env, 1, 2*max_stabilizer_weight);
        qubit_to_variable[q] = x;

        model.add(x);
        all_qubit_vars.add(x);
        // Add constraints.
        model.add( max_of_all >= x ); }
    IloAllDiff uniqueness(env, all_qubit_vars);
    model.add(uniqueness);
    // Setup constraints with other stabilizers.
    for (auto s : gr.get_neighbors(s0)) {
        const auto& other_support = s->support;
        const auto& sch = s->qubit_to_time;
        // We only care about neighboring stabilizers if they
        // already have a schedule.
        if (sch.empty()) continue;

        bool sum_is_nonzero = false;
        IloIntExpr ind_sum(env);
        for (auto q : support) {
            if (std::find(other_support.begin(), other_support.end(), q) == other_support.end()) {
                continue;
            }
            IloIntVar x = qubit_to_variables[q];
            uint t = sch[q];
            model.add( x != t );

            if (s0->qubit_to_pauli[q] != s->qubit_to_pauli[q]) {
                IloBoolVar ineq(env);
                model.add(ineq);
                model.add( ineq == (x >= t) );
                ind_sum += ineq;
                sum_is_nonzero = true;
            }
        }
        if (sum_is_nonzero) {
            IloIntVar sum_ind_div_two(env, 0, IloIntMax);
            model.add(sum_ind_div_two);
            model.add( 2*sum_ind_div_two == ind_sum );
        }
    }
    return model;
}


}   // protean
}   // qontra

#endif  // PROTEAN_SCHEDULER_h
