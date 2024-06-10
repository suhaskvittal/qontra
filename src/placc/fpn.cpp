/*
 *  author: Suhas Vittal
 *  date:   29 May 2024
 * */

#include "placc/fpn.h"

#include <vtils/bijective_map.h>
#include <vtils/two_level_map.h>
#include <vtils/utility.h>

#include <deque>
#include <initializer_list>
#include <map>

namespace placc {

using namespace qontra;
using namespace graph;

FPN::FPN(TannerGraph* tgr, int color_to_remove)
    :Graph(),
    data_qubits(),
    inplace_parity_qubits(),
    removed_parity_qubits(),
    parity_qubits(),
    flag_qubits(),
    obs_list(),
    idctr(0)
{
    vtils::BijectiveMap< sptr<tanner::vertex_t>, sptr<fpn_v_t> >
        tanner_fpn_map;
    for (sptr<tanner::vertex_t> tv : tgr->get_vertices()) {
        // Ignore Z checks as their support is identical to X checks. We only need one or the
        // other.
        if (tv->qubit_type == tanner::vertex_t::type::zparity) continue;
        sptr<fpn_v_t> v = make_and_add_vertex(idctr++);
        if (tv->qubit_type == tanner::vertex_t::type::data) {
            v->qubit_type = fpn_v_t::type::data;
            data_qubits.push_back(v);
        } else {
            v->qubit_type = fpn_v_t::type::parity;
            parity_qubits.push_back(v);
        }
        tanner_fpn_map.put(v, tv);
    }
    // Color the checks of the color code.
    std::map<sptr<tanner::vertex_t>, int> color_map;
    int max_color = tgr->compute_check_color_map(color_map);
    if (max_color != 2) {
        std::cerr << "Coloring is not a 3-coloring." << std::endl;
        exit(1);
    }
    // Now add edges.
    for (sptr<fpn_v_t> v : parity_qubits) {
        sptr<tanner::vertex_t> tv = tanner_fpn_map.at(v);
        std::vector<sptr<fpn_v_t>> support;
        for (sptr<tanner::vertex_t> tw : tgr->get_neighbors(tv)) {
            sptr<fpn_v_t> w = tanner_fpn_map.at(tw);
            support.push_back(w);
        }
        int c = color_map.at(tv);

        v->support = support;
        if (c == color_to_remove) {
            // By removed, we mean removed from the architecture. It still exists in the
            // FPN graph as a placeholder.
            removed_parity_qubits.push_back(v);
        } else {
            inplace_parity_qubits.push_back(v);
            for (sptr<fpn_v_t> w : support) {
                make_and_add_edge(v, w);
            }
        }
    }
    // Set up the observables.
    for (auto obs : tgr->get_obs(false)) {
        std::vector<sptr<fpn_v_t>> _obs;
        for (sptr<tanner::vertex_t> tv : obs) {
            _obs.push_back(tanner_fpn_map.at(tv));
        }
        obs_list.push_back(_obs);
    }
}

void
FPN::place_flags() {
    // By the structure of the color codes, we can actually greedily pair up qubits sharing
    // the same inplace plaquettes, as these qubits will share an edge in the color code
    // lattice.
    vtils::TwoLevelMap<sptr<fpn_v_t>, sptr<fpn_v_t>, size_t> inc_matrix;
    for (sptr<fpn_v_t> p : inplace_parity_qubits) {
        for (size_t i = 0; i < p->support.size(); i++) {
            sptr<fpn_v_t> v = p->support.at(i);
            for (size_t j = i+1; j < p->support.size(); j++) {
                sptr<fpn_v_t> w = p->support.at(j);
                inc_matrix[v][w]++;
                inc_matrix[w][v]++;
            }
        }
    }
    std::set<sptr<fpn_v_t>> used;
    for (sptr<fpn_v_t> v : data_qubits) {
        if (used.count(v)) continue;
        if (!inc_matrix.count(v)) continue; // This is a widowed qubit.
        sptr<fpn_v_t> w = nullptr;
        size_t best_score = 0;
        for (const auto& [u, cnt] : inc_matrix.at(v)) {
            if (!used.count(u) && cnt > best_score) {
                w = u;
                best_score = cnt;
            }
        }
        if (w == nullptr) {
            std::cerr << "[ place_flags ] found qubit without partner" << std::endl;
            exit(1);
        }
        // Make flag qubit.
        sptr<fpn_v_t> f = make_and_add_vertex(idctr++);
        f->qubit_type = fpn_v_t::type::flag;
        for (sptr<fpn_v_t> ch : get_common_neighbors({v,w})) {
            ch->flag_usage_map[std::make_pair(v,w)] = f;
            ch->flag_usage_map[std::make_pair(w,v)] = f;
            for (sptr<fpn_v_t> x : {v,w}) {
                delete_edge(get_edge(x,ch));
            }
            make_and_add_edge(ch,f);
        }
        for (sptr<fpn_v_t> x : {v,w}) {
            make_and_add_edge(x,f);
        }
        flag_qubits.push_back(f);
        vtils::insert_all(used, {v,w});
    }
}

void
FPN::place_widowed_qubits() {
    const int CRITERIA = 4; // TODO: set to 4 once everything is settled.

    std::vector<sptr<fpn_v_t>> widowed;
    for (sptr<fpn_v_t> v : get_vertices()) {
        if (get_degree(v) == 0 && v->qubit_type == fpn_v_t::type::data) {
            widowed.push_back(v);
            v->is_widowed = true;
        }
    }
    // Find appropriate ancilla to place them on.
    for (sptr<fpn_v_t> v : get_vertices()) {
        if (widowed.empty()) break;
        if (v->qubit_type == fpn_v_t::type::data) continue;
        if (get_degree(v) < CRITERIA && get_degree(v) > 0) {
            make_and_add_edge(widowed.back(), v);
            widowed.pop_back();
        }
    }
}

void
FPN::compute_cnot_order() {
    // Basic idea: perform a BFS from an initial parity qubit. As we don't need to
    // consider commutativity rules for the color codes (each check is measured
    // separately on each plaquette), we will attach a timestep for each set of data/flag
    // qubits between two inplace plaquettes. Two parity qubits cannot share the same
    // timestep for the same data/flag qubit.
    struct timestep_grp_t {
        std::set<sptr<fpn_v_t>> adjacent_checks;
        std::vector<sptr<fpn_v_t>> elements;
        std::set<size_t> unavailable_times;
    };
    timestep_map.clear();
    // The timestep_map data structure that is a field of the class is different from the 
    // below data structure, which is merely used to track occupied timesteps.
    std::map<sptr<fpn_v_t>, std::set<size_t>> occupied_timestep_map;

    std::deque<sptr<fpn_v_t>> bfs{ inplace_parity_qubits[0] };
    std::set<sptr<fpn_v_t>> visited;

    max_timestep = 0;
    while (bfs.size()) {
        sptr<fpn_v_t> p = bfs.front(); bfs.pop_front();
        if (visited.count(p)) continue;

        std::set<sptr<fpn_v_t>> adj_checks; // Where to traverse next.
        // Make timestep groups.
        std::vector<timestep_grp_t> timestep_grps;
        for (sptr<fpn_v_t> v : get_neighbors(p)) {
            std::set<sptr<fpn_v_t>> _adj;
            for (sptr<fpn_v_t> w : get_neighbors(v)) {
                if (w->qubit_type == fpn_v_t::type::parity) _adj.insert(w);
            }

            bool found = false;
            for (timestep_grp_t& grp : timestep_grps) {
                if (_adj == grp.adjacent_checks) {
                    grp.elements.push_back(v);
                    vtils::insert_range(
                            grp.unavailable_times, occupied_timestep_map[v]);
                    found = true;
                    break;
                }
            }
            // Make new timestep grp if nothing was found.
            if (!found) {
                timestep_grp_t grp;
                grp.adjacent_checks = _adj;
                grp.unavailable_times = occupied_timestep_map[v];
                timestep_grps.push_back(grp);
            }
            // Update adj_checks.
            vtils::insert_range(adj_checks, _adj);
        }
        // Now we need to schedule the timestep groups.
        const size_t ntg = timestep_grps.size();
        for (size_t t = 0; t < 2*ntg && timestep_grps.size() > 0; t++) {
            for (auto it = timestep_grps.begin(); it != timestep_grps.end(); it++) {
                timestep_grp_t& grp = *it;
                if (!grp.unavailable_times.count(t)) {
                    // Schedule the entire group for timestep t.
                    for (sptr<fpn_v_t> v : grp.elements) {
                        occupied_timestep_map[v].insert(t);
                        timestep_map[v][t] = p;
                    }
                    max_timestep = std::max(t, max_timestep);
                    it = timestep_grps.erase(it);
                    break;
                }
            }
        }
        // Traverse to the next set of checks.
        for (sptr<fpn_v_t> q : adj_checks) {
            bfs.push_back(q);
        }
        visited.insert(p);
    }
}

}   // placc
