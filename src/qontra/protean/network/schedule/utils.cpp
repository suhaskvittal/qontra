/*
 *  author: Suhas Vittal
 *  date:   9 January 2024
 * */

#include "qontra/protean/network.h"

#include <vtils/two_level_map.h>

#include <algorithm>

namespace qontra {
namespace protean {

using namespace net;
using namespace graph;
using namespace vtils;

sptr<raw_vertex_t>
Scheduler::test_and_get_other_endpoint_if_ready(
        sptr<raw_vertex_t> endpoint,
        std::vector<sptr<raw_vertex_t>> path) 
{
    for (size_t i = 1; i < path.size(); i++) {
        sptr<raw_vertex_t> rx = path[i-1],
                            ry = path[i];
        sptr<raw_edge_t> re = net_p->raw_connection_network.get_edge(rx, ry);
        if (rx == endpoint || ry == endpoint) {
            if (!is_good_for_current_cycle(rx) || !is_good_for_current_cycle(ry)) return nullptr;
            if (visited_edge_map[re] < BODY_STAGE) {
                // Then the edge is not done.
                return rx == endpoint ? ry : rx;
            }
        } else {
            if (visited_edge_map[re] < BODY_STAGE) {
                return nullptr;
            }
        }
    }
    return nullptr;
}

Scheduler::cx_t
Scheduler::get_next_edge_between(sptr<raw_vertex_t> src, sptr<raw_vertex_t> dst, bool for_x_check, stage_t s) {
    if (for_x_check) return get_next_edge_between(dst, src, false, s);

    cx_return_status = CX_RET_GOOD;
    if (has_contention(src) || has_contention(dst)) {
        return ret_null_and_set_status(CX_RET_CONTENTION);
    }
    // Check if src and dst are directly connected, or are connected via proxy.
    RawNetwork& raw_net = net_p->raw_connection_network;
    if (raw_net.contains(src, dst)) {
        sptr<raw_edge_t> re = raw_net.get_edge(src, dst);
        if (visited_edge_map[re] < s) {
            if (!is_good_for_current_cycle(src) || !is_good_for_current_cycle(dst)) {
                return ret_null_and_set_status(CX_RET_TOO_EARLY);
            }
            return std::make_tuple(src, dst, re);
        }
    } else {
        // Otherwise, they are connected via proxy.
        std::vector<sptr<raw_vertex_t>> path = net_p->get_proxy_walk_path(src, dst);
        for (size_t i = 0; i < path.size(); i++) {
            sptr<raw_edge_t> re = raw_net.get_edge(path[i-1], path[i]);
            if (visited_edge_map[re] < s) {
                if (!is_good_for_current_cycle(path[i-1]) || !is_good_for_current_cycle(path[i])) {
                    return ret_null_and_set_status(CX_RET_TOO_EARLY);
                }
                // Then, return this edge -- it is not done.
                return std::make_tuple(path[i-1], path[i], re);
            }
        }
    }
    // No edge to be done:
    return ret_null_and_set_status(CX_RET_FINISHED);
}

}   // protean
}   // qontra
