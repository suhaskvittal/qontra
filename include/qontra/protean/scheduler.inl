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
    return pv->cycle_role_map.at(rv) <= cycle
            && test_and_set_physical_qubit(rv);
}

inline bool
Scheduler::has_contention(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    return cx_in_use_set.count(pv);
}

inline bool
Scheduler::test_and_set_physical_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    if (active_role_map[pv] == rv)       { return true; }
    if (active_role_map[pv] == nullptr)  { active_role_map[pv] = rv; return true; }
    return false;
}

inline void
Scheduler::release_physical_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    active_role_map[pv] = nullptr;
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
}

inline void
Scheduler::push_back_measurement(std::vector<uint64_t>& qes_operands, sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    qes_operands.push_back(pv->id);
    meas_ctr_map[rv] = meas_ctr++;
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
