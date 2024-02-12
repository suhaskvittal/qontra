/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#ifndef PROTEAN_SCHEDULER_h
#define PROTEAN_SCHEDULER_h

#include "qontra/protean/network.h"

#include <qontra/ext/qes.h>

#include <vtils/bijective_map.h>
#include <vtils/two_level_map.h>

#include <deque>
#include <set>
#include <tuple>

namespace qontra {
namespace protean {

// The Scheduler class is dedicated to building the syndrome extraction
// schedule.
class Scheduler {
public:
    Scheduler(PhysicalNetwork*);

    qes::Program<>  run(void);

    void    build_preparation(qes::Program<>&);
    void    build_body(qes::Program<>&);
    void    build_teardown(qes::Program<>&);
    void    build_measurement(qes::Program<>&);

    size_t  get_measurement_ctr(void);
    size_t  get_measurement_time(sptr<net::raw_vertex_t>);
    std::map<sptr<net::raw_vertex_t>, size_t>   get_meas_ctr_map(void);
private:
    enum class stage_t {
        needs_undo = -2,
        invalid = -1,
        preparation = 0,
        body = 1,
        teardown = 2,
        measurement = 3,
        done = 4
    };

    enum class cx_ret_t {
        ok = 0,
        too_early = 1,
        done = 2,
        contention = 3,
        proxy_occupied = 4
    };

    // tuple: control, target, corresponding edge, will_need_undo, should_be_considered_in_partial_support
    typedef std::tuple<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_edge_t>, bool, bool>
        cx_t;

    // Indirection functions (inlined):
    //
    // Just a level of indirection to avoid the syntax of net_p->raw_connection_network.get_support(...) (e.g,)
    RawNetwork::parity_support_t&           get_support(sptr<net::raw_vertex_t>);
    std::vector<sptr<net::raw_vertex_t>>&   get_proxy_walk_path(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    stage_t&    get_visited_edge_map(sptr<net::raw_edge_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);

    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>>
        compute_schedule(std::set<sptr<net::raw_vertex_t>>& checks_this_stage);

    bool is_proxy_usable(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t> src, sptr<net::raw_vertex_t> dst);
    bool has_contention(sptr<net::raw_vertex_t>);
    bool test_and_set_qubit(sptr<net::raw_vertex_t>);
    bool test_and_set_path(const std::vector<sptr<net::raw_vertex_t>>&);
    bool test_and_set_proxy_ownership(const std::vector<sptr<net::raw_vertex_t>>&);
    void release_qubit(sptr<net::raw_vertex_t>, bool pop_cycle=true);
    void release_path(const std::vector<sptr<net::raw_vertex_t>>&);
    void release_proxy_ownership(const std::vector<sptr<net::raw_vertex_t>>&);

    cx_t ret_null_and_set_status(cx_ret_t);

    void push_back_cx(
            std::vector<uint64_t>&, cx_t, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, bool, stage_t);
    void push_back_measurement(std::vector<uint64_t>&, sptr<net::raw_vertex_t>);

    void update_proxy_info(sptr<net::raw_vertex_t> proxy);
    void perform_proxy_resets(qes::Program<>&);

    // Retrieves the subset of checks whose stage_map is set to the input.
    std::set<sptr<net::raw_vertex_t>> get_checks_at_stage(stage_t);

    cx_t get_next_edge_between(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, bool for_x_check, stage_t);

    // For the given endpoint, this function (1) traverses the path, stopping early if any edge has not been
    // completed, (2) upon arriving at an edge with the endpoint, returns the other endpoint if the edge
    // is not done, (3) returns nullptr if stopped early or no edge satisfying (2) is found.
    sptr<net::raw_vertex_t> test_and_get_other_endpoint_if_ready(
                                    sptr<net::raw_vertex_t> endpoint,
                                    std::vector<sptr<net::raw_vertex_t>> path,
                                    cx_t&);

    // Note that preparation and teardown are the exact same operations, but teardown is reversed.
    //
    // Gets the H operands (which are X-parity or Z-flags) for a check operator.
    std::vector<uint64_t>   prep_tear_get_h_operands(sptr<net::raw_vertex_t>, stage_t);
    // Gets the CX operands to prepare flags for a check operator.
    std::vector<uint64_t>   prep_tear_get_cx_operands(sptr<net::raw_vertex_t>, stage_t);

    // Gets all data qubits that have a CNOT for this check operator.
    std::set<sptr<net::raw_vertex_t>>
        body_get_partial_data_support(sptr<net::raw_vertex_t>);
    // Gets the prev and post data qubit sets given partial_support.
    std::set<sptr<net::raw_vertex_t>>
        get_prev(sptr<net::raw_vertex_t>);
    std::set<sptr<net::raw_vertex_t>>
        get_post(sptr<net::raw_vertex_t>, const std::set<sptr<net::raw_vertex_t>>&);
    // Gets the CX operands according to the data_cnot_schedules.
    std::vector<uint64_t> 
        body_get_cx_operands(sptr<net::raw_vertex_t>, const std::vector<sptr<net::raw_vertex_t>>&, size_t k);

    // Basic data structures: just holds read-only information.
    std::set<sptr<net::raw_vertex_t>>   all_checks;
    // Tracking structures: used to track where each check is when constructing the schedule.
    std::map<sptr<net::raw_vertex_t>, stage_t>  stage_map;
    // visited_edge_map keys are a bit complex: (edge, src, dst). If the edge is
    // not a proxy edge, then the key is just ((x, y), x, y). If the edge is a
    // proxy edge, then src and dst are the endpoints of the proxy path.
    std::map<std::tuple<sptr<net::raw_edge_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>, stage_t>
        visited_edge_map;
    std::map<sptr<net::raw_vertex_t>, stage_t>  h_gate_stage_map;
    std::map<sptr<net::raw_vertex_t>, stage_t>  flag_stage_map;

    // This map marks if a proxy has a nonzero state by indicating which CX path it is currently
    // involved in. If the value of a proxy key is valid, then the proxy is not usable.
    //
    // 0 = reserved for use, 1 = used once, and 2 = used twice and can be freed.
    std::map<sptr<net::raw_vertex_t>, std::tuple<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, int>>
        proxy_occupied_map;

    std::map<sptr<net::phys_vertex_t>, sptr<net::raw_vertex_t>>
        active_role_map;
    std::map<sptr<net::phys_vertex_t>, std::deque<size_t>>
        cycle_role_order_map;

    std::map<sptr<net::raw_vertex_t>, std::set<sptr<net::raw_vertex_t>>> scheduled_data_qubit_map;
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, stage_t>
        data_stage_map;

    // After probing and checking which data qubits need CNOTs, the CNOT that is required is
    // stored here for easy access.
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, cx_t>
        data_cx_map;
    // Tracks which qubits have been used for CXs.
    std::set<sptr<net::phys_vertex_t>>          cx_in_use_set;
    // Tracks the order of measurements.
    std::map<sptr<net::raw_vertex_t>, size_t>   meas_ctr_map;
    size_t                                      meas_ctr;

    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, double>
        running_ind_sum_map;

    size_t cycle;
    cx_ret_t cx_return_status;

    PhysicalNetwork* net_p;
};

}   // protean
}   // qontra

#include "scheduler.inl"

#endif  // PROTEAN_SCHEDULER_h
