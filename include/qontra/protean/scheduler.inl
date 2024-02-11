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

inline Scheduler::stage_t&
Scheduler::get_visited_edge_map(
        sptr<net::raw_edge_t> re,
        sptr<net::raw_vertex_t> src,
        sptr<net::raw_vertex_t> dst) 
{
    auto k = std::make_tuple(re, src, dst);
    if (!visited_edge_map.count(k)) {
        visited_edge_map[k] = stage_t::invalid;
    }
    return visited_edge_map[k];
}

inline bool
Scheduler::is_proxy_usable(
        sptr<net::raw_vertex_t> rprx,
        sptr<net::raw_vertex_t> src,
        sptr<net::raw_vertex_t> dst) 
{
    if (!proxy_occupied_map.count(rprx)) return true;
    auto& entry = proxy_occupied_map[rprx];
    return std::get<0>(entry) == src && std::get<1>(entry) == dst;
}

inline bool
Scheduler::has_contention(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    return cx_in_use_set.count(pv);
}

inline bool
Scheduler::test_and_set_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    // Check if rv is allowed to test the qubit.
    if (pv->cycle_role_map.at(rv) > cycle) {
        return false;
    }
    // Check if the physical qubit is already taken by rv:
    if (active_role_map[pv] == rv && active_role_map[pv] != nullptr) {
        return true;
    }
    // Check if the physical qubit is unlocked, and also check if rv can take it.
    if (active_role_map[pv] == nullptr
        && (cycle_role_order_map[pv].empty()
            || cycle_role_order_map[pv].front() == pv->cycle_role_map.at(rv)))
    {
        active_role_map[pv] = rv;
        return true;
    }
    return false;
}

inline bool
Scheduler::test_and_set_path(const std::vector<sptr<net::raw_vertex_t>>& path) {
    size_t i = 1;
    for ( ; i < path.size()-1; i++) {
        if (!test_and_set_qubit(path.at(i))) goto failed_test_and_set_path;
    }
    return true;
failed_test_and_set_path:
    // If we fail, we need to unset the prior qubits in the path.
    --i;
    for ( ; i >= 1; i--) {
        release_qubit(path.at(i), false);
    }
    return false;
}

inline bool
Scheduler::test_and_set_proxy_ownership(const std::vector<sptr<net::raw_vertex_t>>& path) {
    sptr<net::raw_vertex_t> src = path[0],
                            dst = path.back();
    size_t i = 1;
    for ( ; i < path.size()-1; i++) {
        sptr<net::raw_vertex_t> r = path.at(i);
        if (is_proxy_usable(r, src, dst)) {
            if (!proxy_occupied_map.count(r)) {
                proxy_occupied_map[r] = std::make_tuple(src, dst, 0);
            }
        } else {
            goto failed_test_and_set_proxy_ownership;
        }
    }
    return true;
failed_test_and_set_proxy_ownership:
    --i;
    for ( ; i >= 1; i--) {
        proxy_occupied_map.erase(path.at(i));
    }
    return false;
}

inline void
Scheduler::release_qubit(sptr<net::raw_vertex_t> rv, bool pop_cycle) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    if (active_role_map[pv] == rv) {
        active_role_map[pv] = nullptr;
        if (pop_cycle && !cycle_role_order_map[pv].empty()) {
            cycle_role_order_map[pv].pop_front();
        }
    }
}

inline void
Scheduler::release_path(const std::vector<sptr<net::raw_vertex_t>>& path) {
    for (size_t i = 1; i < path.size()-1; i++) release_qubit(path.at(i));
}

inline void
Scheduler::release_proxy_ownership(const std::vector<sptr<net::raw_vertex_t>>& path) {
    sptr<net::raw_vertex_t> src = path[0],
                            dst = path.back();
    for (size_t i = 1; i < path.size()-1; i++) {
        sptr<net::raw_vertex_t> r = path.at(i);
        auto& entry = proxy_occupied_map[r];
        if (std::get<0>(entry) == src && std::get<1>(entry) == dst) {
            proxy_occupied_map.erase(r);
        }
    }
}

inline Scheduler::cx_t
Scheduler::ret_null_and_set_status(cx_ret_t s) {
    cx_return_status = s;
    return std::make_tuple(nullptr, nullptr, nullptr, false);
}

inline void
Scheduler::push_back_cx(
        std::vector<uint64_t>& qes_operands,
        cx_t cx,
        sptr<net::raw_vertex_t> src,
        sptr<net::raw_vertex_t> dst,
        bool is_x_check,
        stage_t s) 
{
    if (is_x_check) return push_back_cx(qes_operands, cx, dst, src, false, s);

    sptr<net::raw_vertex_t> rx = std::get<0>(cx),
                        ry = std::get<1>(cx);
    sptr<net::raw_edge_t> re = std::get<2>(cx);
    sptr<net::phys_vertex_t> px = net_p->role_to_phys[rx];
    sptr<net::phys_vertex_t> py = net_p->role_to_phys[ry];

    get_visited_edge_map(re, src, dst) = 
        (net_p->raw_connection_network.contains(src, dst) ^ std::get<3>(cx)) ? s : stage_t::needs_undo;

    vtils::push_back_all(qes_operands, {px->id, py->id});
    vtils::insert_all(cx_in_use_set, {px, py});
    update_proxy_info(rx);
    update_proxy_info(ry);
}

inline void
Scheduler::push_back_measurement(std::vector<uint64_t>& qes_operands, sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = net_p->role_to_phys[rv];
    qes_operands.push_back(pv->id);
    meas_ctr_map[rv] = meas_ctr++;
}

inline void
Scheduler::update_proxy_info(sptr<net::raw_vertex_t> prx) {
    if (proxy_occupied_map.count(prx)) {
        std::get<2>(proxy_occupied_map[prx])++;
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
