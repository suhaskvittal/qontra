/* author: Suhas Vittal
 *  date:   31 May 2023
 * */

#include "protean/compiler.h"

namespace qontra {
namespace protean {

#define PRINT_V(id)  (id >> 30) << "|" << ((id >> 24) & 0x3f) << "|" << (id & 0x00ff'ffff)

Compiler::result_t
Compiler::run(const TannerGraph& tanner_graph, bool verbose) {
    // Set state variables to initial values.
    verbosity = verbose;
    compile_round = 0;

    result_t best_result;
    ir_t ir;
    ir.curr_spec = TannerGraph(tanner_graph);

__place:
    if (verbosity) {
        std::cout << "[ place ] ---------------------\n";
    }
    place(ir);
__reduce:
    if (verbosity) {
        std::cout << "\tcost = " << objective(ir.arch) << "\n";
        std::cout << "\tnumber of qubits = " << ir.arch.get_vertices().size() << "\n";
//      print_connectivity(ir.arch);

        std::cout << "[ reduce ] ---------------------\n";
    }
    reduce(ir);
    if (verbosity) {
        std::cout << "\tcost = " << objective(ir.arch) << "\n";
        std::cout << "\tnumber of qubits = " << ir.arch.get_vertices().size() << "\n";
//      print_connectivity(ir.arch);

        std::cout << "[ schedule ] ---------------------\n";
    }
    schedule(ir);
    if (verbosity) {
        std::cout << "\t#ops = " << ir.schedule.size() << "\n";
        print_schedule(ir.schedule);

        std::cout << "[ score ] ---------------------\n";
    }
    score(ir);
    if (verbosity) {
        std::cout << "\tcost = " << ir.score << ", valid = " << ir.valid << "\n";
    }
    // Update result.
    if (!best_result.valid || (ir.valid && ir.score < best_result.score)) {
        best_result.score = ir.score;
        best_result.valid = ir.valid;
        best_result.arch = ir.arch;
        best_result.schedule = ir.schedule;
    } else if (ir.score >= best_result.score) {
        return best_result;
    }
    // Reset ir (except Tanner graph).
    ir.arch = Processor3D();
    ir.schedule.clear();
    ir.score = 0.0;
    ir.valid = false;

    ir.role_to_qubit.clear();
    ir.qubit_to_roles.clear();
    // Perform transformations on Tanner graph.
    compile_round++;
    if (induce(ir)) {
        goto __place;
    } else if (sparsen(ir)) {
        goto __place;
    } else {
        linearize(ir);
        goto __reduce;
    }
}

void
Compiler::place(ir_t& curr_ir) {
    // We design an architecture according to the following rules:
    //  (1) If X has predecessors, then we connect it to other checks
    //      and data qubits as follows:
    //          (a) We first connect it to other preceding parity checks.
    //          (b) We then connect it to gauge qubits.
    //          (c) Finally, we consider data qubits.
    //      such that the data qubits involved in the connected parity and
    //      gauge qubits are nonintersecting.
    //  (2) If X has no predecessors, we directly connect it to its data qubits.
    
    // First create the qubits for each vertex.
    TannerGraph& tanner_graph = curr_ir.curr_spec;
    Processor3D& arch = curr_ir.arch;
    for (auto v : tanner_graph.get_vertices()) {
        proc3d::vertex_t* w = new proc3d::vertex_t;
        w->id = v->id;
        w->processor_layer = 0;
        arch.add_vertex(w);
    }
    // Now make the connections.
    for (auto tv : tanner_graph.get_vertices()) {
        if (tv->qubit_type == tanner::vertex_t::DATA)    continue;
        auto pv = arch.get_vertex(tv->id);

        auto tv_adj = tanner_graph.get_neighbors(tv);
        auto pred = tanner_graph.get_predecessors(tv);
        if (pred.empty() || (compile_round == 0)) {
            // Just connect to the data qubits.
            for (auto tw : tv_adj) {
                auto pw = arch.get_vertex(tw->id);
                proc3d::edge_t* e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch.add_edge(e);
            }
        } else {
            // Scan for other parity checks.            
            std::set<tanner::vertex_t*> unsat_adj(tv_adj.begin(), tv_adj.end());
            for (auto it = pred.begin(); it != pred.end(); ) {
                tanner::vertex_t* tw = *it;
                auto tw_adj = tanner_graph.get_neighbors(tw);
                bool is_subset_of_unsat = is_subset_of(tw_adj, unsat_adj);
                if (!is_subset_of_unsat) {
                    // No reason to connect to this qubit.
                    it = pred.erase(it);
                    continue;
                }
                if (tw->qubit_type & 0x2) { 
                    // It is another parity qubit 
                    std::set<tanner::vertex_t*> unsat_diff;
                    std::set_difference(unsat_adj.begin(), unsat_adj.end(),
                                        tw_adj.begin(), tw_adj.end(),
                                        std::inserter(unsat_diff, unsat_diff.begin()));
                    unsat_adj = unsat_diff;
                    
                    auto pw = arch.get_vertex(tw->id);
                    auto e = new proc3d::edge_t;
                    e->src = (void*)pv;
                    e->dst = (void*)pw;
                    arch.add_edge(e);
                    it = pred.erase(it);
                } else {
                    it++;
                }
            }
            // Go through once more and try to connect to gauge qubits.
            for (auto it = pred.begin(); it != pred.end(); ) {
                tanner::vertex_t* tw = *it;
                auto tw_adj = tanner_graph.get_neighbors(tw);
                bool is_subset_of_unsat = is_subset_of(tw_adj, unsat_adj);
                if (!is_subset_of_unsat) {
                    // No reason to connect to this qubit.
                    it = pred.erase(it);
                    continue;
                }
                // Anything remaining must be a gauge qubit.
                std::set<tanner::vertex_t*> unsat_diff;
                std::set_difference(unsat_adj.begin(), unsat_adj.end(),
                                    tw_adj.begin(), tw_adj.end(),
                                    std::inserter(unsat_diff, unsat_diff.begin()));
                unsat_adj = unsat_diff;
                auto pw = arch.get_vertex(tw->id);
                auto e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch.add_edge(e);
                it = pred.erase(it);
            }
            // If any data qubits remain in the unsat set, connect to them as well.
            for (auto tw : unsat_adj) {
                auto pw = arch.get_vertex(tw->id);
                auto e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch.add_edge(e);
            }
        }
    }
    // Update IR.
    for (auto tv : tanner_graph.get_vertices()) {
        if (tv->qubit_type == tanner::vertex_t::DATA)   continue;   // No need to track data
                                                                    // qubits, as they never
                                                                    // change or gain roles.
        auto pv = arch.get_vertex(tv->id);
        curr_ir.qubit_to_roles[pv] = std::vector<tanner::vertex_t*>();
        curr_ir.qubit_to_roles[pv].push_back(tv);
        curr_ir.role_to_qubit[tv] = pv;
    }
}

void
Compiler::reduce(ir_t& curr_ir) {
    // Modify the architecture by contracting degree-2 non-data qubits and removing
    // zero-degree qubits (these are obviously gauge qubits).
    TannerGraph& tanner_graph = curr_ir.curr_spec;
    Processor3D& arch = curr_ir.arch;

    for (auto pv : arch.get_vertices()) {
        if (!curr_ir.qubit_to_roles.count(pv))    continue;   // This is a data qubit.
        uint deg = arch.get_degree(pv);
        // First check if it is zero-degree.
        if (deg == 0) {
            arch.delete_vertex(pv);
            for (auto tv : curr_ir.qubit_to_roles[pv]) {
                curr_ir.role_to_qubit.erase(tv);
            }
            curr_ir.qubit_to_roles.erase(pv);
            continue;
        }
        // Otherwise, just skip if deg > 2.
        if (deg > 2)    continue;
        auto pv_adj = arch.get_neighbors(pv);
        // When contracting, we must shift the role of the contracted
        // qubit to a non-data qubit. Check if it has any adjacent
        // non-data qubits.
        //
        // Also check if this qubit is only connected to data qubits. If it is a
        // gauge qubit, then remove it.
        proc3d::vertex_t* victim = nullptr;
        bool any_nondata_neighbors = false;
        for (auto pw : pv_adj) {
            auto tw = tanner_graph.get_vertex(pw->id);
            if (tw->qubit_type != tanner::vertex_t::DATA) {
                victim = pw;
                any_nondata_neighbors = true;
                break;
            }
        }
        if (!any_nondata_neighbors) {
            for (auto tv : curr_ir.qubit_to_roles[pv]) {
                if (tv->qubit_type != tanner::vertex_t::GAUGE)  goto reduce_do_not_remove_non_data;
            }
            // Otherwise delete the qubit.
            arch.delete_vertex(pv);
            for (auto tv : curr_ir.qubit_to_roles[pv]) {
                curr_ir.role_to_qubit.erase(tv);
            }
            curr_ir.qubit_to_roles.erase(pv);
            continue;
        }
reduce_do_not_remove_non_data:
        if (victim == nullptr)  continue;
        if (verbosity) {
            std::cout << "\tReduced vertex " << PRINT_V(pv->id) 
                    << " to " << PRINT_V(victim->id) << "\n";
        }
        for (auto tv : curr_ir.qubit_to_roles[pv]) {
            curr_ir.qubit_to_roles[victim].push_back(tv);
            curr_ir.role_to_qubit[tv] = victim;
        }
        // Delete pv and connect its neighbors to the victim.
        std::vector<proc3d::edge_t*> new_edges; // We're temporarily storing the new edges in a
                                                // std::vector because adding edges will cause
                                                // planarity checks. To not do these prematurely,
                                                // we will allocate the edges and add them after
                                                // deleting pv.
        for (auto pw : pv_adj) {
            if (pw == victim)   continue;
            auto e1 = arch.get_edge(pv, pw);
            auto e2 = new proc3d::edge_t;
            e2->src = (void*)victim;
            e2->dst = (void*)pw;
            new_edges.push_back(e2);
        }
        arch.delete_vertex(pv);
        curr_ir.qubit_to_roles.erase(pv);
        for (auto e : new_edges) {
            arch.add_edge(e);
        }
    }
}

void
Compiler::schedule(ir_t& curr_ir) {
    TannerGraph& tanner_graph = curr_ir.curr_spec;
    Processor3D& arch = curr_ir.arch;
    // Execute the larger checks before the smaller checks.
    std::vector<tanner::vertex_t*>  all_checks;
    for (auto tv : tanner_graph.get_vertices_by_type(tanner::vertex_t::XPARITY)) {
        all_checks.push_back(tv);
    }
    for (auto tv : tanner_graph.get_vertices_by_type(tanner::vertex_t::ZPARITY)) {
        all_checks.push_back(tv);
    }

    typedef struct {
        bool operator()(tanner::vertex_t* x, tanner::vertex_t* y) {
            return g.get_degree(x) > g.get_degree(y);
        }

        TannerGraph& g;
    } check_cmp;
    check_cmp c = {tanner_graph};
    std::sort(all_checks.begin(), all_checks.end(), c);
    // We will use Dijkstra's algorithm to get a table of shortest paths from each data qubit
    // to any check. When trying to execute a check, we will attempt to merge paths from
    // multiple data qubits to compute the CNOT schedule.
    using namespace graph;

    dijkstra::ewf_t<proc3d::vertex_t> idw = [&] (proc3d::vertex_t* v, proc3d::vertex_t* w)
    {
        return 1;
    };
    
    typedef struct {
        std::vector<proc3d::vertex_t*> path;
        bool valid;
    } entry_t;

    dijkstra::callback_t<proc3d::vertex_t, entry_t> d_cb = 
    [&] (proc3d::vertex_t* v1, proc3d::vertex_t* v2,
        const std::map<proc3d::vertex_t*, fp_t>& dist,
        const std::map<proc3d::vertex_t*, proc3d::vertex_t*>& pred)
    {
        entry_t x;
        if (v2 == pred.at(v2)) {
            x.valid = false;
            return x;
        }
        
        auto curr = v2;
        while (curr != v1) {
            x.path.push_back(curr);
            curr = pred.at(curr);
        }
        x.path.push_back(v1);
        std::reverse(x.path.begin(), x.path.end());
        x.valid = true;
        return x;
    };

    std::vector<qc::Instruction> schedule;
    std::deque<qc::Instruction> schedule_buffer;
    
    // To schedule the CNOTs for measuring each parity check, we will create graphs
    // (schedule_graph_t) and perform BFS on them. The idea is that these graphs should be
    // acycle such that the sinks are data qubits and the source is the parity qubit being
    // measured. Then, we only need to reverse the visiting order to figure out the CNOT
    // schedule.

    typedef Graph<proc3d::vertex_t, base::edge_t> schedule_graph_t;

    // set_x_parity and curr_parity_check are state variables that are used to adopt
    // the callback to the graph without making a callback for every parity check.
    bool set_x_parity = false;
    proc3d::vertex_t* curr_parity_check;
    search::callback_t<proc3d::vertex_t> s_cb =
    [&] (proc3d::vertex_t* v1, proc3d::vertex_t* v2)
    {
        qc::Instruction inst;
        inst.name = "CX";
        if (set_x_parity)   inst.operands = std::vector<uint>{v1->id, v2->id};
        else                inst.operands = std::vector<uint>{v2->id, v1->id};
        schedule_buffer.push_front(inst);

        // We also need to include undo operations at the end.
        if (v1 != curr_parity_check)    schedule_buffer.push_back(inst);
    };

    auto distance_matrix = create_distance_matrix(&arch, idw, d_cb);
    for (auto tv : all_checks) {
        schedule_graph_t scheduling_graph;
        std::set<proc3d::vertex_t*> non_data_qubits;

        auto pv = curr_ir.role_to_qubit[tv];
        auto tv_adj = tanner_graph.get_neighbors(tv);
        for (auto td : tv_adj) {
            auto pd = arch.get_vertex(td->id);
            if (!distance_matrix[pd][pv].valid) {
                std::cerr << "ERROR: cannot schedule check " << PRINT_V(tv->id)
                        << " [@" << PRINT_V(pv->id) << "] because data qubit " 
                        << PRINT_V(td->id) << " [@" << PRINT_V(pd->id) << "] has no path "
                        << "to the parity check.\n";
                continue;
            }
            auto path = distance_matrix[pd][pv].path;
            scheduling_graph.add_vertex(pd);
            for (uint i = 1; i < path.size(); i++) {
                auto px = path[i-1];
                auto py = path[i];

                scheduling_graph.add_vertex(py);
                auto e = new base::edge_t;
                e->src = (void*)py; // Reverse the order for BFS
                e->dst = (void*)px;
                scheduling_graph.add_edge(e, false);    // And make sure the edge is directed.
                non_data_qubits.insert(py);
            }
        }
        // Set callback state variables.
        set_x_parity = tv->qubit_type == tanner::vertex_t::XPARITY;
        curr_parity_check = pv;
        // Create H and Mrc gates beforehand
        qc::Instruction h;
        qc::Instruction meas;
        qc::Instruction reset;

        h.name = std::string("H");
        for (auto p : non_data_qubits) {
            h.operands.push_back(p->id);
        }

        meas.name = std::string("Mrc");
        meas.operands.push_back(pv->id);

        reset.name = std::string("R");
        reset.operands.push_back(pv->id);

        // Add the gates to the schedule.
        if (set_x_parity) {
            schedule.push_back(h);
        }
        for (auto e : scheduling_graph.get_edges()) {
            auto src = (tanner::vertex_t*)e->src;
            auto dst = (tanner::vertex_t*)e->dst;
        }
        xfs(&scheduling_graph, pv, s_cb, false);
        for (auto inst : schedule_buffer) schedule.push_back(inst);
        if (set_x_parity) {
            schedule.push_back(h);
        }
        schedule.push_back(meas);
        schedule.push_back(reset);
        // Clear the scheduling buffer for the next qubit.
        schedule_buffer.clear();
    }
    curr_ir.schedule = schedule;
}

void
Compiler::score(ir_t& curr_ir) {
    curr_ir.valid = true;
    for (auto constraint : constraints) {
        curr_ir.valid &= constraint(curr_ir.arch);
    }
    if (curr_ir.valid)  curr_ir.score = objective(curr_ir.arch);
}

bool
Compiler::induce(ir_t& curr_ir) {
    TannerGraph& tanner_graph = curr_ir.curr_spec;
    auto vertices = tanner_graph.get_vertices();

    bool any_change = false;
    
    uint new_max_cw = max_induced_check_weight;
    for (uint i = 0; i < vertices.size(); i++) {
        auto tv = vertices[i];
        if (tv->qubit_type == tanner::vertex_t::DATA
            || tanner_graph.get_degree(tv) > max_induced_check_weight)  continue;
        for (uint j = i+1; j < vertices.size(); j++) {
            auto tw = vertices[j];
            if (tw->qubit_type == tanner::vertex_t::DATA
                || tanner_graph.get_degree(tw) > max_induced_check_weight)  continue;
            auto ti = tanner_graph.induce_predecessor(tv, tw);
            if (ti != nullptr) {
                uint d = tanner_graph.get_degree(ti);
                if (d < new_max_cw) new_max_cw = d;
                any_change = true;
            }
        }
    }
    max_induced_check_weight = new_max_cw;

    return any_change;
}

bool
Compiler::sparsen(ir_t& curr_ir) {
    static uint TANNER_GAUGE_INDEX = 0;

    TannerGraph& tanner_graph = curr_ir.curr_spec;
    auto vertices = tanner_graph.get_vertices();

    bool any_change = false;
    for (auto tv : vertices) {
        if (tv->qubit_type == tanner::vertex_t::DATA)    continue;
        auto tv_adj = tanner_graph.get_neighbors(tv);
        for (uint i = 1; i < tv_adj.size()-1; i++) {
            // We leave up to a degree of 3 (as leaving a degree of 2 will simply
            // result in reduction in the "reduce" pass).
            auto tx = tv_adj[i-1];
            auto ty = tv_adj[i];
            // Create a new gauge qubit.
            tanner::vertex_t* tg = new tanner::vertex_t;
            tg->id = (TANNER_GAUGE_INDEX++) | (tanner::vertex_t::GAUGE << 30);
            tg->qubit_type = tanner::vertex_t::GAUGE;
            tanner_graph.add_vertex(tg);
            // Add corresponding edges.
            tanner::edge_t* tx_tg = new tanner::edge_t;
            tx_tg->src = (void*)tx;
            tx_tg->dst = (void*)tg;
            tanner_graph.add_edge(tx_tg);

            tanner::edge_t* ty_tg = new tanner::edge_t;
            ty_tg->src = (void*)ty;
            ty_tg->dst = (void*)tg;
            tanner_graph.add_edge(ty_tg);

            tanner::edge_t* tg_tv = new tanner::edge_t;
            tg_tv->src = (void*)tg;
            tg_tv->dst = (void*)tv;
            tanner_graph.add_edge(tg_tv);
            // Remove old edges.
            tanner_graph.delete_edge(tanner_graph.get_edge(tx, tg));
            tanner_graph.delete_edge(tanner_graph.get_edge(ty, tg));

            any_change = true;
        }
    }

    return any_change;
}

void
Compiler::linearize(ir_t& curr_ir) {
}

void
print_connectivity(Processor3D& arch) {
    std::cout << "Connections:\n";
    for (auto v : arch.get_vertices()) {
        std::cout << "\tQubit " << PRINT_V(v->id) << ":\t";
        for (auto w : arch.get_neighbors(v)) {
            std::cout << " " << PRINT_V(w->id);
            if (w->is_tsv_junction())   std::cout << "(V)";
        }
        std::cout << "\n";
    }
}

void
print_schedule(const std::vector<qc::Instruction>& sched) {
    std::cout << "Schedule:\n";
    for (auto op : sched) {
        std::cout << "\t" << op.name;
        for (uint x : op.operands) {
            std::cout << " " << PRINT_V(x);
        }
        std::cout << "\n";
    }
}

}   // protean
}   // qontra
