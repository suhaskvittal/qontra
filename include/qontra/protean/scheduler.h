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



    RawNetwork::parity_support_t&           get_support(sptr<net::raw_vertex_t>);

    size_t cycle;

    PhysicalNetwork* net_p;
};

}   // protean
}   // qontra

#include "scheduler.inl"

#endif  // PROTEAN_SCHEDULER_h
