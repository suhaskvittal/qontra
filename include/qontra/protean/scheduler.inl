/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#include <vtils/utility.h>

namespace qontra {
namespace protean {

inline void
Scheduler::push_back_measurement(std::vector<int64_t>& operands, sptr<net::raw_vertex_t> rv) {
    meas_ctr_map[rv] = mctr++;
    operands.push_back(qu(rv));
}

inline void
Scheduler::declare_event_for_qubit(sptr<net::raw_vertex_t> rv) {
    int64_t m = meas_ctr_map.at(rv);
    bool is_cross_round = rv->is_check();
    bool is_memory_x = (rv->qubit_type == net::raw_vertex_t::type::xparity)
                        || (rv->qubit_type == net::raw_vertex_t::type::flag 
                                && !raw_network->x_flag_set.count(rv));
    event_queue.push_back({rv, {m}, is_cross_round, is_memory_x});
}

template <class FUNC> bool
Scheduler::try_and_push_back_cx_operands(
        std::vector<int64_t>& cx_operands,
        std::set<int64_t>& in_use,
        const std::vector<sptr<net::raw_vertex_t>>& path,
        size_t k,
        FUNC additional_test,
        bool reverse_cx_dir)
{
    if (k >= path.size()) return false;
    sptr<net::raw_vertex_t> rx = path.at(k-1),
                            ry = path.at(k);
    if (!additional_test(rx, ry)) return false;
    int64_t qx = qu(rx),
            qy = qu(ry);
    // Check if there is a physical edge between qx and qy.
    sptr<net::phys_vertex_t> px = network->get_vertex(qx),
                             py = network->get_vertex(qy);
    if (!network->contains(px, py)) {
        std::cerr << "[ try_and_push_back_cx_operands ] attempted to do CX(" << qx << ", " << qy << ") when"
            << " no coupling exists." << std::endl;
        exit(1);
    }
    if (in_use.count(qx) || in_use.count(qy)) return false;
    if (reverse_cx_dir) std::swap(qx, qy);
    
    vtils::push_back_all(cx_operands, {qx, qy});
    vtils::insert_all(in_use, {qx, qy});

    return true;
}

inline RawNetwork::parity_support_t&
Scheduler::get_support(sptr<net::raw_vertex_t> v) {
    return network->raw_connection_network->get_support(v);
}

inline int64_t
Scheduler::qu(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = network->role_to_phys.at(rv);
    return static_cast<int64_t>(pv->id);
}

inline bool
Scheduler::test_and_set_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = network->role_to_phys.at(rv);
    if (!active_role_map.count(pv) || active_role_map.at(pv) == nullptr) {
        active_role_map[pv] = rv;
        return true;
    } else if (active_role_map.at(pv) == rv) {
        return true;
    } else {
        return false;
    }
}

inline bool
Scheduler::release_qubit(sptr<net::raw_vertex_t> rv) {
    sptr<net::phys_vertex_t> pv = network->role_to_phys.at(rv);
    if (active_role_map.at(pv) == rv) {
        active_role_map[pv] = nullptr;
        return true;
    }
    return false;
}

inline void
Scheduler::test_and_set_exit_on_fail(sptr<net::raw_vertex_t> rv, std::string caller) {
    if (!test_and_set_qubit(rv)) {
        print_test_and_set_debug_and_exit(rv, caller);
    }
}

inline void
Scheduler::print_test_and_set_debug_and_exit(sptr<net::raw_vertex_t> rv, std::string caller) {
    sptr<net::phys_vertex_t> pv = network->role_to_phys.at(rv);
    std::cerr << "[ " << caller << " ] qubit " << graph::print_v(rv)
        << " failed to acquire " << graph::print_v(pv)
        << " from " << graph::print_v(active_role_map.at(pv))
        << std::endl;
    exit(1);
}

}   // protean
}   // qontra
