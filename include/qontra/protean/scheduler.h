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

    void build_preparation(qes::Program<>&);
    void build_body(qes::Program<>&);
    void build_teardown(qes::Program<>&);

    size_t  get_measurement_ctr(void);
    size_t  get_measurement_time(sptr<net::raw_vertex_t>);
    std::map<sptr<net::raw_vertex_t>, size_t>   get_meas_ctr_map(void);
private:
    // source, target, is_for_x_check
    typedef std::tuple<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>> cx_t;

    void prep_tear_h_gates(qes::Program<>&);
    void prep_tear_cx_gates(qes::Program<>&);

    void schedule_cx_along_path(const std::vector<cx_t>&, qes::Program<>&);
    void push_back_measurement(std::vector<int64_t>&, sptr<net::raw_vertex_t>);

    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>>
        compute_schedules(void);

    template <class FUNC>
    bool try_and_push_back_cx_operands(
            std::vector<int64_t>& cx_operands,
            std::set<int64_t>& in_use,
            const std::vector<sptr<net::raw_vertex_t>>& path,
            size_t k,
            FUNC additional_test,
            bool reverse_cx_dir=false);

    RawNetwork::parity_support_t&   get_support(sptr<net::raw_vertex_t>);
    int64_t                         qu(sptr<net::raw_vertex_t>);

    bool test_and_set_qubit(sptr<net::raw_vertex_t>);
    bool release_qubit(sptr<net::raw_vertex_t>);

    void test_and_set_exit_on_fail(sptr<net::raw_vertex_t>, std::string caller);
    void print_test_and_set_debug_and_exit(sptr<net::raw_vertex_t>, std::string caller);

    std::map<sptr<net::phys_vertex_t>, sptr<net::raw_vertex_t>>
        active_role_map;
    std::vector<sptr<net::raw_vertex_t>>
        checks_this_cycle;

    std::map<sptr<net::raw_vertex_t>, size_t>
        meas_ctr_map;
        
    size_t cycle;
    size_t mctr;

    graph::TannerGraph* tanner_graph;
    uptr<RawNetwork>& raw_network;
    PhysicalNetwork* network;
};

}   // protean
}   // qontra

#include "scheduler.inl"

#endif  // PROTEAN_SCHEDULER_h
