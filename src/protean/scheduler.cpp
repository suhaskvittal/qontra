/* author: Suhas Vittal date:   23 October 2023
 * */

#include "protean/scheduler.h"

#define PROTEAN_OBS_OPT

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
        v->qubit_to_pauli = q2p; gr.add_vertex(v);

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
    // Initialize some of the css_code_data_t structure here.
    css_code_data_t code_data;
    // We will maintain an assignment of tanner graph vertices to qubits.
    std::map<tanner::vertex_t*, uint64_t> tanner_vertex_to_qubit;
    uint64_t next_qubit_id = 0;
    for (auto dv : tanner_graph.get_vertices_by_type(tanner::vertex_t::Type::data)) {
        code_data.data_qubits.push_back(next_qubit_id);
        tanner_vertex_to_qubit[dv] = next_qubit_id;
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
        LPManager<uint>* mgr = construct_scheduling_program(s, gr, &code_data, max_stab_weight);
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
    // Now set the parity qubit information.
    code_data.schedule_depth = max_time;
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

css_code_data_t
make_fault_tolerant(css_code_data_t code_data) {
    bool _g_use_mpi = experiments::G_USE_MPI;
    experiments::G_USE_MPI = false;

    std::vector<stabilizer_t> stab = get_stabilizers(code_data);
    // We will insert a flag for each pair of data qubits in the schedule.
    // i.e. (q1, q2) will get a flag, (q3, q4) will get a flag.
    uint flag_qubit_ctr = code_data.data_qubits.size() + code_data.parity_qubits.size();
    auto add_flags = [&] (std::vector<uint>& parity_list, std::vector<uint>& flag_list)
    {
        for (uint pq : parity_list) {
            const auto& data_qubit_order = code_data.check_schedules[pq];
            int32_t prev = -1;
            for (int32_t dq : data_qubit_order) {
                if (dq < 0) continue;
                if (prev >= 0) {
                    // Then we make a flag.
                    uint fq = flag_qubit_ctr++;
                    code_data.flag_qubits.push_back(fq);
                    flag_list.push_back(fq);
                    if (!code_data.parity_to_flags.count(pq)) {
                        code_data.parity_to_flags[pq] = std::vector<uint>();
                    }
                    code_data.parity_to_flags[pq].push_back(fq);
                    tlm::put(code_data.flag_usage, pq, (uint)prev, fq);
                    tlm::put(code_data.flag_usage, pq, (uint)dq, fq);
                    prev = -1;
                } else {
                    prev = dq;
                }
            }
        }
    };
    add_flags(code_data.zparity_list, code_data.zflag_list);
    add_flags(code_data.xparity_list, code_data.xflag_list);

    experiments::G_USE_MPI = _g_use_mpi;
    return code_data;
}

LPManager<uint>*
construct_scheduling_program(
        stab_vertex_t* s0,
        StabilizerGraph& gr,
        css_code_data_t* css_data_p,
        int max_stabilizer_weight)
{
    LPManager<uint>* mgr = new LPManager<uint>;

    // This variable is greater than all qubit variables.
    lp_var_t<uint> max_of_all = mgr->add_slack_var(1, 2*max_stabilizer_weight, VarBounds::both, VarDomain::integer);

    const auto& support = s0->support;
    // Create variables for each qubit.
    for (uint q : support) {
        lp_var_t<uint> qv = mgr->add_var(q, 1, 2*max_stabilizer_weight, VarBounds::both, VarDomain::integer);
        lp_constr_t<uint> con(max_of_all - qv, 0.0, ConstraintDirection::ge);
        mgr->add_constraint(con);
    }
    // Add uniqueness constraint (any qi != qj).
    const uint w = support.size();
    for (uint i = 0; i < w; i++) {
        lp_var_t<uint> q1v = mgr->get_var(support[i]);
        for (uint j = i+1; j < w; j++) {
            lp_var_t<uint> q2v = mgr->get_var(support[j]);
            lp_constr_t<uint> con(q1v, q2v, ConstraintDirection::neq);
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
            lp_var_t<uint> qv = mgr->get_var(q);
            int t = sch.at(q);
            lp_constr_t<uint> con1(qv, t, ConstraintDirection::neq);
            mgr->add_constraint(con1);
            if (s0->qubit_to_pauli[q] != s->qubit_to_pauli[q]) {
                const double M = 100000;
                lp_var_t<uint> ind = mgr->add_slack_var(0, 1, VarBounds::both, VarDomain::binary);
                // x - t <= My
                lp_constr_t<uint> con2(
                        qv - t, M*ind, ConstraintDirection::le);
                // t - x <= M(1-y)
                lp_constr_t<uint> con3(
                        t - qv, M - M*ind, ConstraintDirection::le);
                mgr->add_constraint(con2);
                mgr->add_constraint(con3);
                // Also add ind to ind_sum.
                ind_sum += lp_expr_t<uint>(ind);
                sum_is_nonzero = true;
            }
        }
        if (sum_is_nonzero) {
            lp_var_t<uint> sum_ind_div_two = mgr->add_slack_var(0, 0, VarBounds::lower, VarDomain::integer);
            lp_constr_t<uint> con(ind_sum, 2*sum_ind_div_two, ConstraintDirection::eq);
            mgr->add_constraint(con);
        }
    }
    // The objective function is as follows (maximization):
    //  -max_of_all
    //      +   sum of (|q_i - q_{i+1}|) normalized by number of pairs.   
    //              if observable optimization is enabled (PROTEAN_OBS_OPT)
    // Build the LP.
    lp_expr_t<uint> obj = -max_of_all;
#ifdef PROTEAN_OBS_OPT
    if (css_data_p != nullptr) {
        // This optimization really only works for CSS codes. I'm not sure how to extend it to general
        // stabilizer codes.
        bool is_x_type = s0->stabilizer[0].first == pauli::x;  // A quick and dirty way to figure out what stabilizer we're working with.
        // X stabilizers may propagate X hooks (and Z is the similarly so). So we care about the same type stabilizer.
        auto obs_ref = is_x_type ? css_data_p->x_obs_list : css_data_p->z_obs_list;

        uint pairs = 0;
        lp_expr_t<uint> obs_opt_obj;
        std::set<std::pair<uint, uint>> visited;
        for (auto obs : obs_ref) {
            for (uint q1 : obs) {
                auto q1v = mgr->get_var(q1);
                for (uint q2 : obs) {
                    if (q1 == q2) continue;

                    auto q1_q2 = std::make_pair(q1, q2);
                    auto q2_q1 = std::make_pair(q2, q1);
                    if (visited.count(q1_q2) || visited.count(q2_q1)) continue;

                    auto q2v = mgr->get_var(q2);
                    pairs++;
                    // We need to add a slack variable and a constraint.
                    auto y = mgr->add_slack_var(0, 0, VarBounds::lower, VarDomain::integer);
                    auto b = mgr->add_slack_var(0, 1, VarBounds::both, VarDomain::binary);
                    // Want y to be max(q1 - q2, q2 - q1)
                    const fp_t M = 100000;
                    lp_constr_t<uint> con1(y, q1v, ConstraintDirection::ge);
                    lp_constr_t<uint> con2(y, q2v, ConstraintDirection::ge);
                    lp_constr_t<uint> con3(y, q1v + M*b, ConstraintDirection::le);
                    lp_constr_t<uint> con4(y, q2v + M*(1.0-b), ConstraintDirection::le);

                    mgr->add_constraint(con1);
                    mgr->add_constraint(con2);
                    mgr->add_constraint(con3);
                    mgr->add_constraint(con4);

                    obs_opt_obj += y;
                    pairs++;
                    visited.insert(q1_q2);
                    visited.insert(q2_q1);
                }
            }
        }
        std::cout << "pairs: " << pairs << "\n";
        obs_opt_obj *= 1.0/((fp_t)pairs);

        obj += obs_opt_obj;
    }
#endif
    mgr->build(obj, true);
    return mgr;
}

}   // protean
}   // qontra
