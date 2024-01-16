/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#include <vtils/utility.h>

namespace qontra {
namespace protean {

inline size_t
Scheduler::get_measurement_ctr() {
    return meas_ctr;
}

inline size_t
Scheduler::get_measurement_time(sptr<net::raw_vertex_t> rv) {
    return meas_ctr_map.at(rv);
}

inline std::map<sptr<net::raw_vertex_t>, size_t>
Scheduler::get_meas_ctr_map() {
    return meas_ctr_map;
}

inline RawNetwork::parity_support_t&
Scheduler::get_support(sptr<net::raw_vertex_t> v) {
    return net_p->raw_connection_network.get_support(v);
}

inline std::vector<sptr<net::raw_vertex_t>>&
Scheduler::get_proxy_walk_path(sptr<net::raw_vertex_t> rx, sptr<net::raw_vertex_t> ry) {
    return net_p->raw_connection_network.get_proxy_walk_path(rx, ry);
}

inline bool
Scheduler::is_good_for_current_cycle(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    return pv->cycle_role_map.at(rv) <= cycle && test_and_set_physical_qubit(rv);
}

inline bool
Scheduler::is_proxy_usable(
        sptr<net::raw_vertex_t> rprx,
        sptr<net::raw_vertex_t> src,
        sptr<net::raw_vertex_t> dst) 
{
    if (!proxy_occupied_map.count(rprx)) return true;
    auto& entry = proxy_occupied_map[rprx];
    return (std::get<0>(entry) == src && std::get<1>(entry) == dst) || (std::get<2>(entry) == 0);
}

inline bool
Scheduler::has_contention(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    return cx_in_use_set.count(pv);
}

inline bool
Scheduler::test_and_set_physical_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    // Check if the physical qubit is already taken by rv:
    if (active_role_map[pv] == rv)       { return true; }
    // Check if the physical qubit is unlocked, and also check if rv can take it.
    if (active_role_map[pv] == nullptr && remaining_roles_map[pv].size() && remaining_roles_map[pv].back() == rv) {
        active_role_map[pv] = rv;
        return true;
    }
    return false;
}

inline void
Scheduler::release_physical_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    if (active_role_map[pv] == rv) {
#ifdef DEBUG
        std::cout << "released " << graph::print_v(rv) << "\n";
#endif
        active_role_map[pv] = nullptr;
        remaining_roles_map[pv].pop_back();
    }
}

inline Scheduler::cx_t
Scheduler::ret_null_and_set_status(int s) {
    cx_return_status = s;
    return std::make_tuple(nullptr, nullptr, nullptr);
}

inline void
Scheduler::push_back_cx(std::vector<uint64_t>& qes_operands, cx_t cx, stage_t s) {
    sptr<net::raw_vertex_t> rx = std::get<0>(cx),
                        ry = std::get<1>(cx);
    sptr<net::raw_edge_t> re = std::get<2>(cx);
    sptr<net::phys_vertex_t> px = net_p->role_to_phys[rx];
    sptr<net::phys_vertex_t> py = net_p->role_to_phys[ry];

    visited_edge_map[re] = s;

    vtils::push_back_all(qes_operands, {px->id, py->id});
    vtils::insert_all(cx_in_use_set, {px, py});
    // Update proxy_reset_map
#ifdef PROTEAN_DEBUG
    std::cerr << "committing CX(" << graph::print_v(rx) << ", " << graph::print_v(ry) << ")" << std::endl;
#endif
    if (proxy_reset_map.count(rx))      proxy_reset_map[rx]++;
    if (proxy_reset_map.count(ry))      proxy_reset_map[ry]++;
    if (proxy_occupied_map.count(rx))   std::get<2>(proxy_occupied_map[rx])++;
    if (proxy_occupied_map.count(ry))   std::get<2>(proxy_occupied_map[ry])++;
}

inline void
Scheduler::push_back_measurement(std::vector<uint64_t>& qes_operands, sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    qes_operands.push_back(pv->id);
    meas_ctr_map[rv] = meas_ctr++;
}

inline void
Scheduler::perform_proxy_resets(qes::Program<>& program) {
    std::vector<uint64_t> operands;
    for (auto it = proxy_reset_map.begin(); it != proxy_reset_map.end(); ) {
        if (it->second == 0) {
            it = proxy_reset_map.erase(it);
        } else if (it->second == 1) {
#ifdef PROTEAN_DEBUG
            std::cerr << "failed reset of proxy " << graph::print_v(it->first) << std::endl;
#endif
            it++;
        } else {
            sptr<net::phys_vertex_t> pprx = net_p->role_to_phys[it->first];
            operands.push_back(pprx->id);

            proxy_occupied_map.erase(it->first);
#ifdef PROTEAN_DEBUG
            std::cerr << "finished reset of proxy " << graph::print_v(it->first) << std::endl;
#endif
            it = proxy_reset_map.erase(it);
        }
    }
    if (operands.size()) {
        program.emplace_back("reset", operands);
    }
    // Also remove any entries in proxy_occupied_map that are <2.
    for (auto it = proxy_occupied_map.begin(); it != proxy_occupied_map.end(); ) {
        if (std::get<2>(it->second) < 2) {
            it = proxy_occupied_map.erase(it);
        } else {
            it++;
        }
    }
}

inline std::set<sptr<net::raw_vertex_t>>
Scheduler::get_checks_at_stage(stage_t s) {
    std::set<sptr<net::raw_vertex_t>> stage_checks(all_checks);
    for (auto it = stage_checks.begin(); it != stage_checks.end(); ) {
        if (stage_map[*it] != s)    it = stage_checks.erase(it);
        else                        it++;
    }
    return stage_checks;
}

}   // protean
}   // qontra
