/*
 *  author: Suhas Vittal
 *  date:   2 June 2024
 * */

#include "placc/cx.h"
#include "placc/fpn.h"
#include "placc/tree.h"

#include <qontra/tables.h>

#include <vtils/utility.h>

#include <algorithm>

using namespace qontra;
using namespace graph;

static tables::ErrorAndTiming et;

namespace placc {

fp_t    // Returns latency under a default ErrorAndTiming model.
write_to_sch(
    qes::Program<>& sch,
    const std::deque<sptr<fpn_v_t>>& commit_queue,
    const std::map<sptr<fpn_v_t>, uptr<ShorTree>>& tree_map,
    bool mx)
{
    // Compute operations for each tree.
    CXManager cxm;

    std::vector<sptr<fpn_v_t>> ancilla, 
                                hstart,
                                hend;
    for (const auto& [x, st] : tree_map) {
        for (sptr<fpn_v_t> v : st->get_vertices()) {
            if (v->qubit_type != fpn_v_t::type::data) {
                ancilla.push_back(v);
            }
        }
        sptr<fpn_v_t> cen = st->get_head();
        hstart.push_back(cen);

        std::deque<sptr<fpn_v_t>> bfs{ cen };
        std::map<sptr<fpn_v_t>, sptr<fpn_v_t>> prev_map;
        std::set<sptr<fpn_v_t>> visited{ cen };
        while (bfs.size()) {
            sptr<fpn_v_t> v = bfs.front(); bfs.pop_front();
            bool all_neighbors_are_data = true;
            for (sptr<fpn_v_t> w : st->get_neighbors(v)) {
                if (visited.count(w)) continue;
                if (w->qubit_type == fpn_v_t::type::data) {
                    if (mx) {
                        cxm.push_back_cx(v->id, w->id);
                    } else {
                        cxm.push_back_cx(w->id, v->id);
                    }
                } else {
                    bfs.push_back(w);
                    cxm.push_back_cx(v->id, w->id);
                    all_neighbors_are_data = false;
                }
            }
            if (prev_map.count(v)) {
                sptr<fpn_v_t> p = prev_map.at(v);
                cxm.push_back_cx(v->id, p->id);
            }
            if (all_neighbors_are_data) {
                hend.push_back(v);
            }
            visited.insert(v);
        }
    }
    fp_t lat = 3*et.t_g1q + cxm.get_depth()*et.t_g2q + et.t_ro;

    push_back_gate(sch, "reset", ancilla);
    push_back_gate(sch, "h", hstart);
    cxm.flush(sch);
    push_back_gate(sch, "h", hend);
    push_back_gate(sch, "measure", ancilla);
    return lat;
}


qes::Program<>
FPN::phase_two_schedule(fp_t& lat) {
    qes::Program<> sch;

    std::set<sptr<fpn_v_t>> blocked_qubits;
    std::deque<sptr<fpn_v_t>> commit_queue;
    std::map<sptr<fpn_v_t>, uptr<ShorTree>> shor_tree_map;

    std::vector<sptr<fpn_v_t>> remaining(removed_parity_qubits);
    size_t subp = 0;
    while (remaining.size()) {
        DMAT m = compute_distance_matrix(blocked_qubits);
        std::map<sptr<fpn_v_t>, uptr<ShorTree>> loc_shor_tree_map;
        std::map<sptr<fpn_v_t>, size_t> cost_map;
        for (sptr<fpn_v_t> x : remaining) {
            sptr<fpn_v_t> cen = get_center_of(x->support, m, blocked_qubits);
            if (cen != nullptr) {
                uptr<ShorTree> shor_tree = std::make_unique<ShorTree>(
                    ShorTree::build_tree(this, cen, x->support, blocked_qubits)
                );
                cost_map[x] = shor_tree->get_spacetime_cost();
                loc_shor_tree_map[x] = std::move(shor_tree);
                std::cout << "[ p2 ] could made tree for " << print_v(x) << " [";
                for (sptr<fpn_v_t> u :x->support) std::cout << " " << print_v(u);
                std::cout << " ], ST-cost=" << cost_map[x] << std::endl;
            } else {
                std::cout << "[ p2 ] could not make tree for " << print_v(x) << std::endl;
                loc_shor_tree_map[x] = nullptr;
                cost_map[x] = std::numeric_limits<size_t>::max();
            }
        }
        // Sort remaining checks by the spacetime cost of their trees.
        std::sort(remaining.begin(), remaining.end(),
                [&] (auto x, auto y)
                {
                    return cost_map[x] < cost_map[y];
                });
        // Greedily select the maximum trees that can be done at once.
        std::set<sptr<fpn_v_t>> new_blocked_qubits;
        bool any_new_commits = false;
        for (auto it = remaining.begin(); it != remaining.end(); ) {
            auto x = *it;
            if (loc_shor_tree_map.at(x) == nullptr) { it++; continue; }
            uptr<ShorTree>& st = loc_shor_tree_map.at(x);
            // Check if any newly blocked qubits are in the tree. If so, skip.
            bool found_blockage = false;
            for (sptr<fpn_v_t> v : new_blocked_qubits) {
                if (st->contains(v)) {
                    found_blockage = true;
                    std::cout << "\tfailed to schedule " << print_v(x) << " due to blockage on " << print_v(v) << std::endl;
                    break;
                }
            }
            if (found_blockage) { it++; continue; }
            // Otherwise, this qubit can be committed.
            std::cout << "\t" << print_v(x) << " induced blocks on";
            for (sptr<fpn_v_t> v : st->get_vertices()) std::cout << " " << print_v(v);
            std::cout << std::endl;
            vtils::insert_range(new_blocked_qubits, st->get_vertices());
            shor_tree_map[x] = std::move(st);
            commit_queue.push_back(x);
            it = remaining.erase(it);
            any_new_commits = true;
        }
        std::cout << "[ p2 ] remaining =";
        for (sptr<fpn_v_t> x : remaining) std::cout << " " << print_v(x);
        std::cout << std::endl;
        if (any_new_commits) {
            vtils::insert_range(blocked_qubits, new_blocked_qubits);
            continue;
        }
        // Otherwise, attempt to schedule CNOTs.
        std::cout << "[ p2 ] subphase " << subp << ", cumlat = " << lat << "ns" << std::endl;
        subp++;
        for (int mx = 0; mx <= 1; mx++) {
            lat += write_to_sch(sch, commit_queue, shor_tree_map, static_cast<bool>(mx));
        }
        commit_queue.clear();
        shor_tree_map.clear();
        blocked_qubits.clear();
    }
    // Final commit:
    std::cout << "[ p2 ] subphase " << subp << ", cumlat = " << lat << "ns" << std::endl;
    subp++;
    for (int mx = 0; mx <= 1; mx++) {
        lat += write_to_sch(sch, commit_queue, shor_tree_map, static_cast<bool>(mx));
    }
    return sch;
}

}   // placc
