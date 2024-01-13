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
    typedef int stage_t;
    typedef std::tuple<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_edge_t>>
        cx_t;

    // Indirection functions (inlined):
    //
    // Just a level of indirection to avoid the syntax of net_p->raw_connection_network.get_support(...) (e.g,)
    RawNetwork::parity_support_t&           get_support(sptr<net::raw_vertex_t>);
    std::vector<sptr<net::raw_vertex_t>>&   get_proxy_walk_path(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);

    bool    is_good_for_current_cycle(sptr<net::raw_vertex_t>);
    bool    has_contention(sptr<net::raw_vertex_t>);
    bool    test_and_set_physical_qubit(sptr<net::raw_vertex_t>);
    void    release_physical_qubit(sptr<net::raw_vertex_t>);

    cx_t    ret_null_and_set_status(int);

    void    push_back_cx(std::vector<uint64_t>&, cx_t, stage_t);
    void    push_back_measurement(std::vector<uint64_t>&, sptr<net::raw_vertex_t>);

    // Retrieves the subset of checks whose stage_map is set to the input.
    std::set<sptr<net::raw_vertex_t>>       get_checks_at_stage(stage_t);

    cx_t get_next_edge_between(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, bool for_x_check, stage_t);

    // For the given endpoint, this function (1) traverses the path, stopping early if any edge has not been
    // completed, (2) upon arriving at an edge with the endpoint, returns the other endpoint if the edge
    // is not done, (3) returns nullptr if stopped early or no edge satisfying (2) is found.
    sptr<net::raw_vertex_t> test_and_get_other_endpoint_if_ready(
                                    sptr<net::raw_vertex_t> endpoint, std::vector<sptr<net::raw_vertex_t>> path);

    // Note that preparation and teardown are the exact same operations, but teardown is reversed.
    //
    // Gets the H operands (which are X-parity or Z-flags) for a check operator.
    std::vector<uint64_t>   prep_tear_get_h_operands(sptr<net::raw_vertex_t>, stage_t);
    // Gets the CX operands to prepare flags for a check operator.
    std::vector<uint64_t>   prep_tear_get_cx_operands(sptr<net::raw_vertex_t>, stage_t);

    // Gets all data qubits that have a CNOT for this check operator.
    std::set<sptr<net::raw_vertex_t>>
        body_get_partial_data_support(sptr<net::raw_vertex_t>);
    // Gets the CX operands according to the data_cnot_schedules.
    std::vector<uint64_t> 
        body_get_cx_operands(sptr<net::raw_vertex_t>, const std::vector<sptr<net::raw_vertex_t>>&, size_t k);

    // Basic data structures: just holds read-only information.
    std::set<sptr<net::raw_vertex_t>>   all_checks;
    // Tracking structures: used to track where each check is when constructing the schedule.
    std::map<sptr<net::raw_vertex_t>, stage_t>  stage_map;
    std::map<sptr<net::raw_edge_t>, stage_t>    visited_edge_map;
    std::map<sptr<net::raw_vertex_t>, stage_t>  h_gate_stage_map;
    std::map<sptr<net::raw_vertex_t>, stage_t>  flag_stage_map;

    std::map<sptr<net::phys_vertex_t>, sptr<net::raw_vertex_t>>
        active_role_map;

    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, stage_t>
        data_stage_map;

    // After probing and checking which data qubits need CNOTs, the CNOT that is required is
    // stored here for easy access.
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, cx_t>
        data_cx_map;
    // Tracks which qubits have been used for CXs.
    std::set<sptr<net::phys_vertex_t>>      cx_in_use_set;
    // Tracks the order of measurements.
    std::map<sptr<net::raw_vertex_t>, size_t>   meas_ctr_map;
    size_t                                      meas_ctr;

    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, double>
        running_ind_sum_map;
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, size_t>
        running_ind_ctr_map;

    size_t  cycle;
    int     cx_return_status;

    PhysicalNetwork* net_p;

    const int CX_RET_GOOD = 0;
    const int CX_RET_TOO_EARLY = 1;
    const int CX_RET_FINISHED = 2;
    const int CX_RET_CONTENTION = 3;

    const stage_t PREP_STAGE = 0;
    const stage_t BODY_STAGE = 1;
    const stage_t TEAR_STAGE = 2;
    const stage_t MEAS_STAGE = 3;
    const stage_t DONE_STAGE = 4;
};

}   // protean
}   // qontra

#include "scheduler.inl"

#endif  // PROTEAN_SCHEDULER_h
