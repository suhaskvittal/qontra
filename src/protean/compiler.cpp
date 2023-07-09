/* author: Suhas Vittal 
 * date:   31 May 2023 
 * */

#include "protean/compiler.h"

namespace qontra {
namespace protean {

static uint64_t MIDDLEMAN_INDEX = 1;

using namespace graph;
using namespace compiler;

ir_t*
Compiler::run(TannerGraph* tanner_graph, std::string sdl_file) {
    // Set state variables to initial values.
    compile_round = 0;

    uint rounds_without_progress = 0;

    // Initialize IR.
    SDGraph ideal_sch = build_schedule_graph_from_sdl(filename);
    ideal_sch->dealloc_on_delete = false;

    ir_t* best_result = nullptr;
    ir_t* ir = new ir_t;
    ir->curr_spec = tanner_graph;
    ir->arch = new Processor3D;
    ir->schedule_graph = SDGraph(ideal_sch);
__place:
    if (params.verbose) {
        std::cout << "[ place ] ---------------------\n";
    }
    place(ir);
    if (params.verbose) {
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_mean_connectivity() << ", max = "
                    << ir->arch->get_max_connectivity() << "\n";
        std::cout << "\tthickness = " << ir->arch->get_thickness() << "\n";
        std::cout << "[ unify ] ---------------------\n";
    }
    unify(ir);
__reduce:
    if (params.verbose) {
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_mean_connectivity() << ", max = "
                    << ir->arch->get_max_connectivity() << "\n";
        std::cout << "\tthickness = " << ir->arch->get_thickness() << "\n";
        std::cout << "[ reduce ] ---------------------\n";
    }
    reduce(ir);
    // Check if we need to call merge or split
    if (params.verbose) {
        std::cout << "\tnumber of qubits = " << ir->arch->get_vertices().size() << "\n";
        std::cout << "\tconnectivity = " << ir->arch->get_mean_connectivity() << ", max = "
                    << ir->arch->get_max_connectivity() << "\n";
        std::cout << "\tthickness = " << ir->arch->get_thickness() << "\n";
    }
__constraints:
    if (rounds_without_progress > 0) {
        if (check_size_violation(ir)) {
            if (params.verbose) std::cout << "[ merge ] ---------------------\n";
            if (merge(ir))  goto __reduce;
        }
        if (check_connectivity_violation(ir)) {
            if (params.verbose) std::cout << "[ split ] ---------------------\n";
            if (split(ir))  goto __reduce;
        }
        if (check_thickness_violation(ir)) {
            if (params.verbose) std::cout << "[ flatten ] ------------------------\n";
            if (flatten(ir))    goto __reduce;
        }
    }
__schedule:
    if (params.verbose) {
        std::cout << "[ schedule ] ---------------------\n";
    }
    xform_schedule(ir);
    if (params.verbose) {
        std::cout << "\t#ops = " << ir->schedule.size() << "\n";
//      std::cout << "\tdepth = " << ir->dependency_graph->get_depth() << "\n";

        std::cout << "[ score ] ---------------------\n";
    }
__score:
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
    } else if (ir->score <= best_result->score + 1e-2) {
        std::cout << "\tNO PROGRESS.\n";
        rounds_without_progress++;
        if (rounds_without_progress >= 10) return best_result;
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
    new_ir->schedule_graph = SDGraph(ideal_sch);
    // Perform transformations on Tanner graph.
    compile_round++;
    if (params.verbose)      std::cout << "[ induce ] ---------------------\n";
    if (induce(new_ir)) {
        if (best_result != ir) delete ir;
        ir = new_ir;
        goto __place;
    }

    if (params.verbose)      std::cout << "[ sparsen ] ---------------------\n";
    if (sparsen(new_ir)) {
        if (best_result != ir) delete ir;
        ir = new_ir;
        goto __place;
    }

    delete new_ir;

    if (params.verbose)      std::cout << "[ raise ] -----------------------\n";
    raise(ir);
    goto __constraints;
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
Compiler::unify(ir_t* curr_ir) {
    // No need to check for constraint violations in unify
    // as unify will always reduce the number of qubits and
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
        const auto& checks = pair.second;
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
    // Reallocate edges on the processor.
    curr_ir->arch->reallocate_edges();
}

void
Compiler::reduce(ir_t* curr_ir) {
    // We need not check if any connectivity constraints are violated as connectivity is
    // maintained.
    //
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
        if (curr_ir->is_data(pv))   continue;
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
            if (tw == nullptr)  continue;
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
            if (tw == nullptr)  continue;
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
    arch->reallocate_edges();
}

bool
Compiler::merge(ir_t* curr_ir) {
    Processor3D* arch = curr_ir->arch;

    uint min_total_degree = std::numeric_limits<uint>::max();
    proc3d::edge_t* min_edge = nullptr;
    for (auto e : arch->get_edges()) {
        auto v1 = (proc3d::vertex_t*)e->src;
        auto v2 = (proc3d::vertex_t*)e->dst;
        // Make sure neither qubit is a data qubit.
        if (curr_ir->is_data(v1) || curr_ir->is_data(v2))   continue;
        uint d = arch->get_degree(v1) + arch->get_degree(v2);
        if (d < min_total_degree) {
            min_total_degree = d;
            min_edge = e;
        }
    }
    if (min_edge == nullptr)    return false;
    // Merge the qubits in the min_edge.
    auto v1 = (proc3d::vertex_t*)min_edge->src;
    auto v2 = (proc3d::vertex_t*)min_edge->dst;
    std::vector<proc3d::edge_t*> new_edges;
    for (auto w : arch->get_neighbors(v2)) {
        if (w == v1)    continue;
        if (arch->contains(v1, w))  continue;

        auto f = new proc3d::edge_t;
        f->src = (void*)v1;
        f->dst = (void*)w;
        new_edges.push_back(f);
    }
    // Transfer roles of v2 to v1 and delete v2.
    for (auto tv : curr_ir->qubit_to_roles[v2]) {
        curr_ir->qubit_to_roles[v1].push_back(tv);
        curr_ir->role_to_qubit[tv] = v1;
    }
    if (params.verbose) {
        std::cout << "\tDeleted qubit " << PRINT_V(v2->id) << " and merged it with "
                << PRINT_V(v1->id) << ".\n";
    }
    arch->delete_vertex(v2);
    // Add new edges
    for (auto f : new_edges)    arch->add_edge(f);
    return true;
}

bool
Compiler::split(ir_t* curr_ir) {
    Processor3D* arch = curr_ir->arch;
    
    uint max_degree = 0;
    proc3d::vertex_t* target;
    for (auto v : arch->get_vertices()) {
        uint d = arch->get_degree(v);
        if (d > max_degree) {
            target = v;
            max_degree = d;
        }
    }
    if (max_degree <= 3)    return false; // Don't even try -- no point.

    // Create the partner qubit
    proc3d::vertex_t* dup = new proc3d::vertex_t;
    dup->id = ((target->id & ID_MASK) << ID_GEN_OFFSET) 
                | MIDDLEMAN_INDEX++ 
                | (tanner::vertex_t::GAUGE << ID_TYPE_OFFSET);
    arch->add_vertex(dup);
    auto target_dup = new proc3d::edge_t;
    target_dup->src = target;
    target_dup->dst = dup;
    // Create a dummy entry for roles
    curr_ir->qubit_to_roles[dup] = std::vector<tanner::vertex_t*>();
    // Create adjacency list -- equally partition neighbors of target between target and dup.
    auto tar_adj = arch->get_neighbors(target);

    std::vector<proc3d::edge_t*> new_edges;
    for (uint i = 0; i < (max_degree>>1); i++) {
        auto pv = tar_adj[i];
        proc3d::edge_t* e = new proc3d::edge_t;
        e->src = (void*)dup;
        e->dst = (void*)pv;
        new_edges.push_back(e);
        // Delete the edge between target and pv.
        auto f = arch->get_edge(target, pv);
        arch->delete_edge(f);
    }

    arch->add_edge(target_dup);
    for (auto e : new_edges)    arch->add_edge(e);
    return true;
}

bool
Compiler::flatten(ir_t* curr_ir) {
    Processor3D* arch = curr_ir->arch;
    proc3d::edge_t* violator = nullptr;
    for (auto e : arch->get_edges()) {
        auto impl = arch->get_physical_edges(e);
        if (impl.size() > 1) {
            uint layer = impl[1]->processor_layer;
            if (layer > constraints.max_thickness) {
                violator = e;
                break;
            }
        }
    }
    if (violator == nullptr)    return false;

    proc3d::vertex_t* src = (proc3d::vertex_t*)violator->src;
    proc3d::vertex_t* dst = (proc3d::vertex_t*)violator->dst;
    proc3d::vertex_t* mm = new proc3d::vertex_t;
    mm->id = ((src->id & ID_MASK) << ID_GEN_OFFSET)
                | MIDDLEMAN_INDEX++
                | (tanner::vertex_t::GAUGE << ID_TYPE_OFFSET);
    curr_ir->qubit_to_roles[mm] = std::vector<tanner::vertex_t*>();

    if (params.verbose) {
        std::cout << "\tCreated gauge " << PRINT_V(mm->id) << " for coupling between "
            << PRINT_V(src->id) << " and " << PRINT_V(dst->id) << ".\n";
    }
    
    auto src_mm = new proc3d::edge_t;
    src_mm->src = src;
    src_mm->dst = mm;
    auto mm_dst = new proc3d::edge_t;
    mm_dst->src = mm;
    mm_dst->dst = dst;

    arch->delete_edge(violator);
    arch->add_edge(src_mm);
    arch->add_edge(mm_dst);
    return true;
}

void
Compiler::xform_schedule(ir_t* curr_ir) {
    // Transform the existing schedule graph so it can be run on
    // the defined coupling graph.
    TannerGraph* tanner_graph = curr_ir->curr_spec;
    Processor3D* arch = curr_ir->arch;
    SDGraph& sch_graph = curr_ir->schedule_graph;

    distance::callback_t d_cb = [&] (proc3d::vertex_t* v1, 
                                    proc3d::vertex_t* v2,
                                    const std::map<proc3d::vertex_t*, fp_t> dist,
                                    const std::map<proc3d::vertex_t*, proc3d::vertex_t*> prev)
    {
        std::vector<proc3d::vertex_t*> path;
        auto curr = v2;
        while (curr != v1) {
            path.push_back(curr);
            curr = prev[v2];
        }
        path.push_back(v1);
        return path;
    };
    auto path_table = distance::create_distance_matrix(
                                arch, 
                                unit_ewf_t<proc3d::vertex_t>(),
                                d_cb);

    std::map<proc3d::vertex_t*, uint32_t>   pv_to_qubitno;
    uint32_t pvid = 0;
    for (auto pv : arch->get_vertices()) {
        pv_to_qubitno[pv] = pvid++;
    }

    for (auto sdv : sch_graph.get_vertices()) {
        auto sch = sdv->sch;
        schedule_t new_sch;
        schedule_t epilogue;    // Any extra instructions we may need
                                // to execute
        std::map<proc3d::vertex_t*, uint> gauge_acc_ctr;
        for (auto& inst : sch) {
            if (inst.name == "cx") {
                uint dq = inst.operands[0];
                uint pq = inst.operands[1];
                if ((dq >> 30) & 0x2)   std::swap(dq, pq);
                bool is_x_check = (pq >> 30) == tanner::vertex_t::XPARITY;

                // Note that sdv's have 32-bit ids that are used in the
                // instructions.
                uint64_t pq64 = (pq & ((1<<30)-1)) | ((pq >> 30) << ID_TYPE_OFFSET);
                auto tdv = tanner_graph->get_vertex(dq);
                auto tpv = tanner_graph->get_vertex(pq);
                // Get corresponding physical qubits.
                auto pdv = arch->get_vertex(tdv);
                auto ppv = curr_ir->role_to_qubit[tpv];

                auto path = path_table[pdv][ppv];
                for (uint i = 1; i < path.size(); i++) {
                    auto px = path[i-1];
                    auto py = path[i];
                    if ((py >> ID_TYPE_OFFSET) == tanner::vertex_t::GAUGE) {
                        // Only proceed if all data qubits 
                    }
                }
            }
        }
    }
}


void
Compiler::score(ir_t* curr_ir) {
    curr_ir->arch->reallocate_edges();
    // Check if constraints are maintained.
    curr_ir->valid = !check_connectivity_violation(curr_ir)
                && !check_size_violation(curr_ir)
                && curr_ir->arch->get_thickness() <= constraints.max_thickness;
    // If constraints are maintained, score the IR.
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

    uint max_degree = 0;
    tanner::vertex_t* target = nullptr;
    // Get maximum degree vertex in tanner graph.
    for (auto tv : vertices) {
        if (tv->qubit_type == tanner::vertex_t::DATA
            || curr_ir->sparsen_visited_set.count(tv))  continue;
        uint d = tanner_graph->get_degree(tv);
        if (d > max_degree) {
            target = tv;
            max_degree = d;
        }
    }
    if (target == nullptr)  return false;
    auto target_adj = tanner_graph->get_neighbors(target);
    for (uint i = 1; i < target_adj.size()-1; i += 2) {
        // We leave up to a degree of 3 (as leaving a degree of 2 will simply
        // result in reduction in the "reduce" pass).
        auto tx = target_adj[i-1];
        auto ty = target_adj[i];
        std::vector<tanner::vertex_t*> tg_adj{tx, ty};
        if (tanner_graph->has_copy_in_gauges(tg_adj)) continue;
        // Create a new gauge qubit.
        tanner::vertex_t* tg = new tanner::vertex_t;
        uint64_t index = tanner_graph->get_vertices_by_type(tanner::vertex_t::GAUGE).size();
        tg->id = (index) | (tanner::vertex_t::GAUGE << ID_TYPE_OFFSET);
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
            std::cout << "\t( " << PRINT_V(target->id) << " ) Added new gauge between " 
                            << PRINT_V(tx->id) << " and "
                            << PRINT_V(ty->id) << "\n";
        }
    }
    return true;
}

void
Compiler::raise(ir_t* curr_ir) {
    Processor3D* arch = curr_ir->arch;
    auto cpl_length_table = arch->get_coupling_lengths();
    std::vector<proc3d::edge_t*> in_plane_edges;
    for (auto e : arch->get_edges()) {
        if (!arch->has_complex_coupling(e)) in_plane_edges.push_back(e);
    }

    auto cmp = [&] (proc3d::edge_t* e1, proc3d::edge_t* e2) {
        return cpl_length_table[e1] > cpl_length_table[e2];
    };
    auto longest_coupling = *(std::min_element(
                                    in_plane_edges.begin(),
                                    in_plane_edges.end(),
                                    cmp));
    arch->force_out_of_plane(longest_coupling);
}

stim::Circuit
build_stim_circuit(
        compiler::ir_t* ir, 
        uint rounds, 
        const std::vector<uint>& obs, 
        bool is_memory_x,
        ErrorTable& errors,
        TimeTable& times) 
{
    auto tanner_graph = ir->curr_spec;
    auto arch = ir->arch;
    auto dependency_graph = ir->dependency_graph;
    // First assign data qubits, then assign other qubits.    
    // This will keep the numbering consistent with the spec provided by the user.
    uint k = 0;
    std::map<uint, uint> id_to_num;
    for (auto v : tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA)) {
        id_to_num[v->id] = v->id;   // Note that data qubits have the same ID as
                                    // in the physical architecture.
        if (v->id > k)  k = v->id;
    }
    const uint n_data = k+1;
    for (auto v : arch->get_vertices()) {
        if (ir->is_data(v)) continue;
        id_to_num[v->id] = ++k;
    }
    const uint n_nondata = k - n_data;
    const uint n = k;

    stim::Circuit circuit;
    // Add initialization for memory experiment.
    for (auto v : arch->get_vertices()) {
        uint i = id_to_num[v->id];
        circuit.append_op("R", {i});
#ifndef DISABLE_ERRORS
        circuit.append_op("X_ERROR", {i}, errors.op1q["R"][v->id]); 
#endif
        if (is_memory_x && i < n_data) {
            circuit.append_op("H", {i});
#ifndef DISABLE_ERRORS
            circuit.append_op("DEPOLARIZE1", {i}, errors.op1q["H"][v->id]); 
#endif
        }
    }
    // Add syndrome extraction rounds to circuit. 
    // Note that the schedule is one round.
    fp_t prev_round_length = 0.0;
    uint prev_round_meas = 0;
    for (uint r = 0; r < rounds; r++) {
        circuit.append_op("TICK", {});
        // Apply decoherence on data qubits.
        for (uint i = 0; i < n_data; i++) {
            fp_t e = 1.0 - exp(-prev_round_length / times.t1[i]);
#ifndef DISABLE_ERRORS
            circuit.append_op("DEPOLARIZE1", {i}, e);
#endif
        }
        // Schedule operations by depth in the circuit.
        prev_round_length = 0.0;

        uint this_round_meas = 0;
        std::vector<uint> meas_events;
        for (uint d = 1; d <= dependency_graph->get_depth(); d++) {
            circuit.append_op("TICK", {});

            auto ops = dependency_graph->get_vertices_at_depth(d);
            fp_t layer_length = 0.0;
            for (auto v : ops) {
                Instruction* inst = v->inst_p;
                if (inst->name == "NOP")    continue;

                std::vector<uint> operands = inst->operands;
                // Add operation and add errors depending on which operation it is.
                fp_t max_opt = 0;
                if (inst->name == "CX") {
                    for (uint j = 0; j < operands.size(); j += 2) {
                        uint x1 = operands[j];
                        uint x2 = operands[j+1];
                        auto x1_x2 = std::make_pair(x1, x2);

                        uint i1 = id_to_num[x1];
                        uint i2 = id_to_num[x2];

                        circuit.append_op(inst->name, {i1, i2});
#ifndef DISABLE_ERRORS
                        circuit.append_op("L_TRANSPORT", {i1, i2},
                                errors.op2q_leakage_transport["CX"][x1_x2]);
                        circuit.append_op("L_ERROR", {i1, i2},
                                errors.op2q_leakage_injection["CX"][x1_x2]);
                        circuit.append_op("DEPOLARIZE2", {i1, i2},
                                errors.op2q["CX"][x1_x2]);
#endif
                        // To implement crosstalk, we will just apply depolarizing
                        // errors on nearby qubits.
                        /*
                        auto pv1 = arch->get_vertex(x1);
                        auto pv2 = arch->get_vertex(x2);
                        for (auto pw : arch->get_neighbors(pv1)) {
                            if (pw == pv2)  continue;
                            uint iw = id_to_num[pw->id];
                            circuit.append_op("DEPOLARIZE1", {iw}, 
                                    errors.op2q_crosstalk["CX"][x1_x2]);
                        }
                        for (auto pw : arch->get_neighbors(pv2)) {
                            if (pw == pv1)  continue;
                            uint iw = id_to_num[pw->id];
                            circuit.append_op("DEPOLARIZE1", {iw}, 
                                    errors.op2q_crosstalk["CX"][x1_x2]);
                        }
                        fp_t t = times.op2q["CX"][x1_x2];
                        if (t > max_opt)    max_opt = t;
                        */
                    }
                } else {
                    for (uint x : operands) {
                        uint i = id_to_num[x];
                        if (inst->name == "Mrc" || inst->name == "Mnrc") {
#ifndef DISABLE_ERRORS
                            circuit.append_op("X_ERROR", {i}, errors.op1q["Mrc"][x]);
#endif
                            circuit.append_op("M", {i});

                            // If this measuring the same check as the memory experiment
                            // targets, then mark it as a recorded detection event.
                            if (inst->is_measuring_x_check == is_memory_x) {
                                meas_events.push_back(this_round_meas);
                            }
                            this_round_meas++;
                        } else if (inst->name == "R") {
                            circuit.append_op(inst->name, {i});
#ifndef DISABLE_ERRORS
                            circuit.append_op("X_ERROR", {i}, errors.op1q["R"][x]);
#endif
                        } else {
                            circuit.append_op(inst->name, {i});
#ifndef DISABLE_ERRORS
                            circuit.append_op("DEPOLARIZE1", {i}, errors.op1q["H"][x]);
#endif
                        }
                        fp_t t = times.op1q[inst->name][x];
                        if (t > max_opt)    max_opt = t;
                    }
                }
                if (layer_length < max_opt) layer_length = max_opt;
            }
            prev_round_length += layer_length;
        }
        // Add error detection events.
        if (r == 0) {
            for (uint m : meas_events) {
                circuit.append_op("DETECTOR", 
                        { (this_round_meas - m) | stim::TARGET_RECORD_BIT });
            }
        } else {
            for (uint m : meas_events) {
                circuit.append_op("DETECTOR", 
                    { (this_round_meas - m) | stim::TARGET_RECORD_BIT,
                        (this_round_meas - m + prev_round_meas) | stim::TARGET_RECORD_BIT });
            }
        }
        prev_round_meas = this_round_meas;
    }
    // Finally, measure all the data qubits.
    for (uint i = 0; i < n_data; i++) {
        // Don't have any errors -- makes life easier.
        if (is_memory_x) {
            circuit.append_op("H", {i});
        }
        circuit.append_op("M", {i});
    }
    // Record observable.
    std::vector<uint> obs_inc;
    for (uint i : obs) {
        obs_inc.push_back((n_data - i) | stim::TARGET_RECORD_BIT);
    }
    std::sort(obs_inc.begin(), obs_inc.end());
    circuit.append_op("OBSERVABLE_INCLUDE", obs_inc, 0);
    return circuit;
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
        spec_out << "GXZ"[tv->qubit_type-1] << (tv->id & ID_MASK);
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
    std::map<uint64_t, uint64_t> id_to_num;
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
        if (impl.size() > 1) {
            d3d_map_out << ",L" << impl[1]->processor_layer;
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
        label_out << (tv->id & ID_MASK) << "," << id_to_num[pv->id] << "\n";
    }
    // Also write data qubits
    for (auto tv : tanner_graph->get_vertices_by_type(tanner::vertex_t::DATA)) {
        label_out << "D" << (tv->id & ID_MASK) << "," << id_to_num[tv->id] << "\n";
    }
    // Write schedule of operations
    std::filesystem::path sched = folder/"schedule.qasm";
    std::ofstream sched_out(sched);

    uint cbit = 0;
    std::string sched_str;
    for (auto inst : ir->schedule) {
        if (inst.name == "NOP") continue;
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
