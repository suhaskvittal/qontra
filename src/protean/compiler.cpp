/* 
 * author: Suhas Vittal 
 * date:   31 May 2023 
 * */

#include "protean/compiler.h"

namespace qontra {
namespace protean {

#define PRINT_V(id)  (id >> 30) << "|" << ((id >> 24) & 0x3f) << "|" << (id & 0x00ff'ffff)

using namespace compiler;

ir_t*
Compiler::run(TannerGraph* tanner_graph) {
    // Set state variables to initial values.
    compile_round = 0;

    uint rounds_without_progress = 0;

    ir_t* best_result = nullptr;
    ir_t* ir = new ir_t;
    ir->curr_spec = tanner_graph;
    ir->arch = new Processor3D;
__place:
    if (params.verbose) {
        std::cout << "[ place ] ---------------------\n";
    }
    place(ir);
    if (params.verbose) {
        std::cout << "\tcost = " << objective(ir) << "\n";
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_connectivity() << "\n";
        std::cout << "[ merge ] ---------------------\n";
    }
    merge(ir);
__reduce:
    if (params.verbose) {
        std::cout << "\tcost = " << objective(ir) << "\n";
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_connectivity() << "\n";
        std::cout << "[ reduce ] ---------------------\n";
    }
    reduce(ir);
    if (params.verbose) {
        std::cout << "\tcost = " << objective(ir) << "\n";
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_connectivity() << "\n";

        std::cout << "[ schedule ] ---------------------\n";
    }
    schedule(ir);
    if (params.verbose) {
        std::cout << "\t#ops = " << ir->schedule.size() << "\n";

        std::cout << "[ score ] ---------------------\n";
    }
    score(ir);
    if (params.verbose) {
        std::cout << "\tcost = " << ir->score << ", valid = " << ir->valid << "\n";
        std::cout << "\trounds without progress = " << rounds_without_progress << "\n";
    }
    // Update result.
    if (best_result == nullptr 
        || !best_result->valid 
        || (ir->valid && ir->score < best_result->score)) 
    {
        best_result = ir;
    } else if (ir->score == best_result->score) {
        rounds_without_progress++;
        if (rounds_without_progress >= 2) return best_result;
    } else if (ir->score > best_result->score) {
        return best_result;
    }
    if (!ir->valid) {
        rounds_without_progress++;
        if (rounds_without_progress >= 10) return best_result;
    }
    // Reset ir (except Tanner graph).
    ir_t* new_ir = new ir_t;
    new_ir->curr_spec = ir->curr_spec;
    new_ir->arch = new Processor3D;
    // Perform transformations on Tanner graph.
    compile_round++;
    if (params.verbose)      std::cout << "[ induce ] ---------------------\n";
    if (induce(new_ir)) {
        delete ir;
        ir = new_ir;
        goto __place;
    }

    if (params.verbose)      std::cout << "[ sparsen ] ---------------------\n";
    if (sparsen(ir)) {
        delete ir;
        ir = new_ir;
        goto __place;
    }
    // If we cannot induce or sparsen to modify the tanner graph, revert back
    // to the previous IR and try to modify connections there.
    delete new_ir;
    if (params.verbose)      std::cout << "[ linearize ] ---------------------\n";
    //linearize(ir);
    goto __reduce;
}

void
Compiler::place(ir_t* curr_ir) {
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
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    Processor3D* arch = curr_ir->arch;
    for (auto v : tanner_graph->get_vertices()) {
        proc3d::vertex_t* w = new proc3d::vertex_t;
        w->id = v->id;
        w->processor_layer = 0;
        arch->add_vertex(w);
    }
    // Now make the connections.
    for (auto tv : tanner_graph->get_vertices()) {
        if (tv->qubit_type == tanner::vertex_t::DATA)    continue;
        auto pv = arch->get_vertex(tv->id);

        auto tv_adj = tanner_graph->get_neighbors(tv);
        auto pred = tanner_graph->get_predecessors(tv);
        if (pred.empty() || (compile_round == 0)) {
            // Just connect to the data qubits.
            for (auto tw : tv_adj) {
                auto pw = arch->get_vertex(tw->id);
                proc3d::edge_t* e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch->add_edge(e);
            }
        } else {
            // Scan for other parity checks.            
            std::set<tanner::vertex_t*> unsat_adj(tv_adj.begin(), tv_adj.end());
            for (auto it = pred.begin(); it != pred.end(); ) {
                tanner::vertex_t* tw = *it;
                auto tw_adj = tanner_graph->get_neighbors(tw);
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
                    
                    auto pw = arch->get_vertex(tw->id);
                    auto e = new proc3d::edge_t;
                    e->src = (void*)pv;
                    e->dst = (void*)pw;
                    arch->add_edge(e);
                    it = pred.erase(it);
                } else {
                    it++;
                }
            }
            // Go through once more and try to connect to gauge qubits.
            for (auto it = pred.begin(); it != pred.end(); ) {
                tanner::vertex_t* tw = *it;
                auto tw_adj = tanner_graph->get_neighbors(tw);
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
                auto pw = arch->get_vertex(tw->id);
                auto e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch->add_edge(e);
                it = pred.erase(it);
            }
            // If any data qubits remain in the unsat set, connect to them as well.
            for (auto tw : unsat_adj) {
                auto pw = arch->get_vertex(tw->id);
                auto e = new proc3d::edge_t;
                e->src = (void*)pv;
                e->dst = (void*)pw;
                arch->add_edge(e);
            }
        }
    }
    // Update IR.
    for (auto tv : tanner_graph->get_vertices()) {
        if (tv->qubit_type == tanner::vertex_t::DATA)   continue;   // No need to track data
                                                                    // qubits, as they never
                                                                    // change or gain roles.
        auto pv = arch->get_vertex(tv->id);
        curr_ir->qubit_to_roles[pv] = std::vector<tanner::vertex_t*>();
        curr_ir->qubit_to_roles[pv].push_back(tv);
        curr_ir->role_to_qubit[tv] = pv;

        if (tv->qubit_type == tanner::vertex_t::GAUGE)  curr_ir->is_gauge_only.insert(pv);
    }
}

void
Compiler::merge(ir_t* curr_ir) {
    // No need to check for constraint violations in merge
    // as merge will always reduce the number of qubits and
    // connectivity.

    TannerGraph* tanner_graph = curr_ir->curr_spec;
    std::map<uint, std::vector<tanner::vertex_t*>> weight_to_checks;
    for (auto tv : tanner_graph->get_vertices()) {
        if (tv->qubit_type != tanner::vertex_t::XPARITY
            && tv->qubit_type != tanner::vertex_t::ZPARITY) continue;
        uint w = tanner_graph->get_degree(tv);
        if (!weight_to_checks.count(w)) weight_to_checks[w] = std::vector<tanner::vertex_t*>();
        weight_to_checks[w].push_back(tv);
    }
    // Now, check if any checks in each equivalence class have the same adjacency list.
    for (auto pair : weight_to_checks) {
        auto checks = pair.second;
        std::set<tanner::vertex_t*> already_removed;
        for (uint i = 0; i < checks.size(); i++) {
            auto tv = checks[i];
            if (already_removed.count(tv))  continue;

            auto tv_adj = tanner_graph->get_neighbors(tv);
            for (uint j = i+1; j < checks.size(); j++) {
                auto tw = checks[j];
                if (already_removed.count(tw))  continue;

                auto tw_adj = tanner_graph->get_neighbors(tw);
                if (tv_adj.size() != tw_adj.size()) continue;
                std::vector<tanner::vertex_t*> intersect_adj;
                std::set_symmetric_difference(tv_adj.begin(), tv_adj.end(),
                                    tw_adj.begin(), tw_adj.end(),
                                    std::back_inserter(intersect_adj));
                if (intersect_adj.empty()) {
                    // Delete the second check and merge the roles.
                    auto pv = curr_ir->role_to_qubit[tv];
                    auto pw = curr_ir->role_to_qubit[tw];
                    if (params.verbose) {
                        std::cout << "\tMerging " << PRINT_V(pw->id) << " with "
                                << PRINT_V(pv->id) << "\n";
                    }
                    curr_ir->role_to_qubit[tw] = pv;
                    
                    for (auto tu : curr_ir->qubit_to_roles[pw]) {
                        curr_ir->qubit_to_roles[pv].push_back(tu);
                    }
                    curr_ir->arch->delete_vertex(pw);
                    curr_ir->qubit_to_roles.erase(pw);

                    already_removed.insert(tw);
                }
            }
        }
    }
}

void
Compiler::reduce(ir_t* curr_ir) {
    // Modify the architecture by contracting degree-2 non-data qubits and removing
    // zero-degree qubits (these are obviously gauge qubits).
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    Processor3D* arch = curr_ir->arch;
    // Solve this as a maximum matching problem.
    using namespace lemon;
    ListGraph matching_graph;
    ListGraph::EdgeMap<int> graph_weights(matching_graph);

    std::map<proc3d::vertex_t*, ListGraph::Node> proc2lemon;
    std::map<ListGraph::Node, proc3d::vertex_t*> lemon2proc;

    for (auto pv : arch->get_vertices()) {
        if (!curr_ir->qubit_to_roles.count(pv))    continue;   // This is a data qubit.
        uint deg = arch->get_degree(pv);
        // First check if it is zero-degree.
        if (deg == 0) {
            arch->delete_vertex(pv);
            for (auto tv : curr_ir->qubit_to_roles[pv]) {
                curr_ir->role_to_qubit.erase(tv);
            }
            curr_ir->qubit_to_roles.erase(pv);
            continue;
        }
        // Otherwise, just skip if deg > 2.
        if (deg > 2)    continue;
        auto pv_adj = arch->get_neighbors(pv);
        // When contracting, we must shift the role of the contracted
        // qubit to a non-data qubit. Check if it has any adjacent
        // non-data qubits.
        //
        // Also check if this qubit is only connected to data qubits. If it is a
        // gauge qubit, then remove it.
        uint victim_degree = std::numeric_limits<uint>::max();
        proc3d::vertex_t* victim = nullptr;
        bool victim_is_gauge = false;
        bool any_nondata_neighbors = false;
        for (auto pw : pv_adj) {
            auto tw = tanner_graph->get_vertex(pw->id);
            if (tw->qubit_type != tanner::vertex_t::DATA)   any_nondata_neighbors = true;
        }
        if (!any_nondata_neighbors) {
            if (!curr_ir->is_gauge_only.count(pv))   goto reduce_do_not_remove_nondata;
            // Otherwise delete the qubit.
            if (params.verbose) {
                std::cout << "\tDeleting qubit " << PRINT_V(pv->id) << "\n";
            }
            for (auto tv : curr_ir->qubit_to_roles[pv]) {
                curr_ir->role_to_qubit.erase(tv);
            }
            curr_ir->qubit_to_roles.erase(pv);
            curr_ir->is_gauge_only.erase(pv);
            arch->delete_vertex(pv);
            continue;
        }
reduce_do_not_remove_nondata:
        // Create node for pv.
        ListGraph::Node pv_node;
        if (!proc2lemon.count(pv)) {
            pv_node = matching_graph.addNode();
            proc2lemon[pv] = pv_node;
            lemon2proc[pv_node] = pv;
        } else {
            pv_node = proc2lemon[pv];
        }

        for (auto pw : pv_adj) {
            auto tw = tanner_graph->get_vertex(pw->id);
            if (tw->qubit_type == tanner::vertex_t::DATA)   continue;
            uint pw_deg = arch->get_degree(pw);
            if (pw_deg <= 2)    continue;
            ListGraph::Node pw_node;
            if (!proc2lemon.count(pw)) {
                pw_node = matching_graph.addNode();
                proc2lemon[pw] = pw_node;
                lemon2proc[pw_node] = pw;
            } else {
                pw_node = proc2lemon[pw];
            }
            auto e = matching_graph.addEdge(pv_node, pw_node);
            // Prioritize matching to gauge qubits.
            graph_weights[e] = curr_ir->is_gauge_only.count(pw) ? 1 : 1;
        }
    }
    // Now that we have the matching structure completed, solve for the maximum matching.
    MaxWeightedMatching matching(matching_graph, graph_weights);
    matching.run();
    std::set<proc3d::vertex_t*> visited;
    std::set<proc3d::vertex_t*> deallocated;
    std::vector<proc3d::edge_t*> new_edges; // We're temporarily storing the new edges in a
                                            // std::vector because adding edges will cause
                                            // planarity checks. To not do these prematurely,
                                            // we will allocate the edges and add them after
                                            // deleting pv.
    for (auto pair : proc2lemon) {
        auto pv = pair.first;
        auto pv_node = pair.second;
        auto pw_node = matching.mate(pv_node);
        auto pw = lemon2proc[pw_node];

        if (pw == nullptr || pv == nullptr)  continue;
        if (visited.count(pv) || visited.count(pw)) continue;
        
        uint deg_pv = arch->get_degree(pv);
        if (deg_pv > 2) {
            // Swap pv and pw
            auto tmp = pv;
            pv = pw; pw = tmp;
        }
        // Delete pv and connect its neighbors to pw.
        auto pv_adj = arch->get_neighbors(pv);
        for (auto pu : pv_adj) {
            if (pu == pw)   continue;
            auto e = new proc3d::edge_t;
            e->src = (void*)pw;
            e->dst = (void*)pu;
            new_edges.push_back(e);
        }
        curr_ir->is_gauge_only.erase(pw);
        for (auto tv : curr_ir->qubit_to_roles[pv]) {
            curr_ir->role_to_qubit[tv] = pw;
            curr_ir->qubit_to_roles[pw].push_back(tv);
        }
        curr_ir->qubit_to_roles.erase(pv);

        if (params.verbose) {
            std::cout << "\tReducing " << PRINT_V(pv->id) << " to " << PRINT_V(pw->id) << "\n";
        }
        
        visited.insert(pv);
        visited.insert(pw);
        deallocated.insert(pv);
    }
    for (auto pv : deallocated) arch->delete_vertex(pv);
    for (auto e : new_edges)    arch->add_edge(e);
}

void
Compiler::schedule(ir_t* curr_ir) {
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    Processor3D* arch = curr_ir->arch;
    // Execute the larger checks before the smaller checks.
    std::vector<tanner::vertex_t*>  all_checks;
    for (auto tv : tanner_graph->get_vertices_by_type(tanner::vertex_t::XPARITY)) {
        all_checks.push_back(tv);
    }
    for (auto tv : tanner_graph->get_vertices_by_type(tanner::vertex_t::ZPARITY)) {
        all_checks.push_back(tv);
    }

    typedef struct {
        bool operator()(tanner::vertex_t* x, tanner::vertex_t* y) {
            return g->get_degree(x) > g->get_degree(y);
        }

        TannerGraph* g;
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

    schedule_t<qc::Instruction> schedule;
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

    auto distance_matrix = create_distance_matrix(arch, idw, d_cb);
    for (auto tv : all_checks) {
        schedule_graph_t scheduling_graph;
        std::set<proc3d::vertex_t*> non_data_qubits;

        auto pv = curr_ir->role_to_qubit[tv];
        auto tv_adj = tanner_graph->get_neighbors(tv);
        for (auto td : tv_adj) {
            auto pd = arch->get_vertex(td->id);
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
                e->is_undirected = false;
                if (!scheduling_graph.add_edge(e))  delete e;
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
        xfs(&scheduling_graph, pv, s_cb, false);
        for (auto inst : schedule_buffer) schedule.push_back(inst);
        if (set_x_parity) {
            schedule.push_back(h);
        }
        schedule.push_back(meas);
        schedule.push_back(reset);
        // Clear the scheduling buffer for the next qubit.
        schedule_buffer.clear();
        // Deallocate entries in scheduling graph manually: we want to destroy
        // the edges but not the vertices.
        scheduling_graph.dealloc_on_delete = false;
        for (auto e : scheduling_graph.get_edges()) delete e;
    }
    curr_ir->schedule = schedule;
}

void
Compiler::score(ir_t* curr_ir) {
    curr_ir->valid = true;
    for (auto constraint : constraints) {
        curr_ir->valid &= constraint(curr_ir);
    }
    if (curr_ir->valid)  curr_ir->score = objective(curr_ir);
}

bool
Compiler::induce(ir_t* curr_ir) {
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    auto vertices = tanner_graph->get_vertices();

    bool any_change = false;
    
    uint new_max_cw = max_induced_check_weight;
    for (uint i = 0; i < vertices.size(); i++) {
        auto tv = vertices[i];
        if (tv->qubit_type == tanner::vertex_t::DATA
            || tanner_graph->get_degree(tv) > max_induced_check_weight)  continue;
        for (uint j = i+1; j < vertices.size(); j++) {
            auto tw = vertices[j];
            if (tw->qubit_type == tanner::vertex_t::DATA
                || tanner_graph->get_degree(tw) > max_induced_check_weight)  continue;
            auto ti = tanner_graph->induce_predecessor(tv, tw);
            if (ti != nullptr) {
                uint d = tanner_graph->get_degree(ti);
                if (d < new_max_cw) new_max_cw = d;
                any_change = true;
                if (params.verbose) {
                    std::cout << "\tInduced gauge " << PRINT_V(ti->id)
                                << " on " << PRINT_V(tv->id) << " and " << PRINT_V(tw->id)
                                << " succeeded.\n";
                }
            }
        }
    }
    max_induced_check_weight = new_max_cw;

    return any_change;
}

bool
Compiler::sparsen(ir_t* curr_ir) {
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    auto vertices = tanner_graph->get_vertices();

    bool any_change = false;
    for (auto tv : vertices) {
        if (tv->qubit_type == tanner::vertex_t::DATA)    continue;
        auto tv_adj = tanner_graph->get_neighbors(tv);
        for (uint i = 1; i < tv_adj.size()-1; i += 2) {
            // We leave up to a degree of 3 (as leaving a degree of 2 will simply
            // result in reduction in the "reduce" pass).
            auto tx = tv_adj[i-1];
            auto ty = tv_adj[i];
            std::vector<tanner::vertex_t*> tg_adj{tx, ty};
            if (tanner_graph->has_copy_in_gauges(tg_adj)) continue;
            // Create a new gauge qubit.
            tanner::vertex_t* tg = new tanner::vertex_t;
            uint index = tanner_graph->get_vertices_by_type(tanner::vertex_t::GAUGE).size();
            tg->id = (index) | (tanner::vertex_t::GAUGE << 30);
            tg->qubit_type = tanner::vertex_t::GAUGE;
            tanner_graph->add_vertex(tg);
            // Add corresponding edges.
            tanner::edge_t* tx_tg = new tanner::edge_t;
            tx_tg->src = (void*)tx;
            tx_tg->dst = (void*)tg;
            tanner_graph->add_edge(tx_tg);

            tanner::edge_t* ty_tg = new tanner::edge_t;
            ty_tg->src = (void*)ty;
            ty_tg->dst = (void*)tg;
            tanner_graph->add_edge(ty_tg);
            
            if (params.verbose) {
                std::cout << "\t( " << PRINT_V(tv->id) << " ) Added new gauge between " 
                                << PRINT_V(tx->id) << " and "
                                << PRINT_V(ty->id) << "\n";
            }

            any_change = true;
        }
    }
    return any_change;
}

void
Compiler::linearize(ir_t* curr_ir) {
    // Here, we will search for more reduction opportunities by searching for parity/gauge qubits
    // with neighboring degree 1 data qubits, such that if we move these qubits, we can reduce
    // the parity/gauge qubit (as it is now degree 2).
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    Processor3D* arch = curr_ir->arch;

    std::vector<proc3d::edge_t*> edges_to_remove;
    std::vector<proc3d::edge_t*> edges_to_add;

    for (auto tv : tanner_graph->get_vertices()) {
        auto pv = arch->get_vertex(tv->id);
        uint deg = arch->get_degree(pv);
        if (deg > 1)    continue;
        // Check if any neighboring parity qubits match the criteria.
        proc3d::vertex_t* victim1 = nullptr;
        for (auto pw : arch->get_neighbors(pv)) {
            if (!curr_ir->qubit_to_roles.count(pw))  continue;   // This is a data qubit.
            uint deg = arch->get_degree(pw);
            if (deg == 3) {
                victim1 = pw;
                break;
            }
        }
        if (victim1 == nullptr)  continue;
        // Search for a data victim adjacent to pw.
        uint victim2_dg = std::numeric_limits<uint>::max();
        proc3d::vertex_t* victim2 = nullptr;
        for (auto pw : arch->get_neighbors(victim1)) {
            if (curr_ir->qubit_to_roles.count(pw))  continue;    // This is not a data qubit.
            if (pw == pv)   continue;
            deg = arch->get_degree(pw);
            if (deg < victim2_dg) {
                victim2_dg = deg;
                victim2 = pw;
            }
        }
        if (victim2 == nullptr) continue;
        // Now, we delete the edge between pv and victim1 and move pv adjacent to victim2.
        auto e1 = arch->get_edge(pv, victim1);
        arch->delete_edge(e1);

        auto e2 = new proc3d::edge_t;
        e2->src = (void*)pv;
        e2->dst = (void*)victim2;
        arch->add_edge(e2);
    }
}

void
print_connectivity(Processor3D* arch) {
    std::cout << "Connections:\n";
    for (auto v : arch->get_vertices()) {
        std::cout << "\tQubit " << PRINT_V(v->id) << ":\t";
        for (auto w : arch->get_neighbors(v)) {
            std::cout << " " << PRINT_V(w->id);
            if (w->is_tsv_junction())   std::cout << "(V)";
        }
        std::cout << "\n";
    }
}

void
print_schedule(const schedule_t<qc::Instruction>& sched) {
    std::cout << "Schedule:\n";
    for (auto op : sched) {
        std::cout << "\t" << op.name;
        for (uint x : op.operands) {
            std::cout << " " << PRINT_V(x);
        }
        std::cout << "\n";
    }
}

void
write_ir_to_folder(ir_t* ir, std::string folder_name) {
    std::filesystem::path folder(folder_name);
    safe_create_directory(folder);
    std::filesystem::path arch_folder = folder/"arch";
    safe_create_directory(arch_folder);
    // Write spec.txt
    std::filesystem::path spec = folder/"spec.txt";
    std::ofstream spec_out(spec);

    TannerGraph* tanner_graph = ir->curr_spec;
    for (auto tv : tanner_graph->get_vertices()) {
        if (tv->qubit_type == tanner::vertex_t::DATA)   continue;
        spec_out << "GXZ"[tv->qubit_type-1] << (tv->id & 0x00ff'ffff);
        for (auto tw : tanner_graph->get_neighbors(tv)) {
            spec_out << "," << tw->id;
        }
        spec_out << "\n";
    }
    // Write arch/coupling files.
    std::filesystem::path flat_map = arch_folder/"flat_map.txt";
    std::filesystem::path d3d_map = arch_folder/"3d_map.txt";

    std::ofstream flat_map_out(flat_map);
    std::ofstream d3d_map_out(d3d_map);
    
    Processor3D* arch = ir->arch;
    // First, organize the qubits into more appropriate ids.
    std::map<uint, uint> id_to_num;
    uint k = 0;
    for (auto v : arch->get_vertices())  id_to_num[v->id] = k++;
    // Now, write the architecture.
    for (auto e : arch->get_edges()) {
        auto pv = (proc3d::vertex_t*)e->src;
        auto pw = (proc3d::vertex_t*)e->dst;
        // Flat map is pretty straightforward.
        flat_map_out << id_to_num[pv->id] << "," << id_to_num[pw->id] << "\n";
        // Check if the connection is complex.
        auto impl = arch->get_physical_edges(pv, pw);
        d3d_map_out << id_to_num[pv->id];
        for (auto f : impl) {
            auto py = (proc3d::vertex_t*)f->dst;
            d3d_map_out << "," << py->id;
        }
        d3d_map_out << "," << id_to_num[pw->id] << "\n";
    }
    // Write Tanner vertex to physical qubit mapping.
    std::filesystem::path labels = folder/"labels.txt";
    std::ofstream label_out(labels);

    auto& role_to_qubit = ir->role_to_qubit;
    for (auto pair : role_to_qubit) {
        auto tv = pair.first;
        auto pv = pair.second;
        label_out << "DGXZ"[tv->qubit_type];
        label_out << (tv->id & 0x00ff'ffff) << "," << id_to_num[pv->id] << "\n";
    }
    // Also write data qubits
    for (auto tv : tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA)) {
        label_out << "D" << (tv->id & 0x00ff'ffff) << "," << id_to_num[tv->id] << "\n";
    }
    // Write schedule of operations
    std::filesystem::path sched = folder/"schedule.qasm";
    std::ofstream sched_out(sched);

    uint cbit = 0;
    std::string sched_str;
    for (auto inst : ir->schedule) {
        if (inst.name == "Mrc") {
            sched_str += "measure q[" + std::to_string(id_to_num[inst.operands[0]]) 
                + "] -> c[" + std::to_string(cbit++) + "];\n";
        } else {
            std::string opname;
            if (inst.name == "R")   opname = "reset";
            else {
                for (auto x : inst.name)    opname.push_back(std::tolower(x));
            }
            uint ops_per_call = opname == "cx" ? 2 : 1;
            for (uint i = 0; i < inst.operands.size(); i += ops_per_call) {
                sched_str += opname + " q[" + std::to_string(id_to_num[inst.operands[i]]) + "]";
                if (ops_per_call == 2) {
                    sched_str += ",q[" + std::to_string(id_to_num[inst.operands[i+1]]) + "]";
                }
                sched_str += ";\n";
            }
        }
    }
    sched_out << "OPENQASM 2.0;\n"
                << "include \"qelib1.inc\";\n"
                << "qreg q[" << arch->get_vertices().size() << "];\n"
                << "creg c[" << cbit << "];\n"
                << sched_str;
}

}   // protean
}   // qontra
