/*
 *  author: Suhas Vittal
 *  date:   31 May 2024
 * */

#include "placc/fpn.h"

#include <qontra/graph/algorithms/distance.h>

namespace placc {

DMAT
FPN::compute_distance_matrix(const std::set<sptr<fpn_v_t>>& blocked_qubits) {
    DMAT m = create_distance_matrix<fpn_v_t, fpn_e_t, fp_t>(this, 
            [&] (sptr<fpn_e_t> e)
            {
                sptr<fpn_v_t> x = e->get_source<fpn_v_t>(),
                              y = e->get_target<fpn_v_t>();
                if (blocked_qubits.count(x) || blocked_qubits.count(y)) {
                    return 1000;
                } else {
                    return 1;
                }
            },
            [&] (sptr<fpn_v_t> x, sptr<fpn_v_t> y, const auto& dist, const auto& pred)
            {
                if (x != y && (!pred.count(y) || pred.at(y) == nullptr)) {
                    return -1.0;
                } else {
                    return dist.at(y);
                }
            });
    return m;
}

sptr<fpn_v_t>
FPN::get_center_of(
        const std::vector<sptr<fpn_v_t>>& qubits,
        const DMAT& m,
        const std::set<sptr<fpn_v_t>>& blocked_qubits)
{
    // Compute the radius of each vertex.
    std::map<sptr<fpn_v_t>, fp_t> rad_map;
    std::set<sptr<fpn_v_t>> inf_rad; // Mark vertices with infinite radius (disconnected).
    for (sptr<fpn_v_t> v : qubits) {
        for (sptr<fpn_v_t> x : get_vertices()) {
            if (x->qubit_type == fpn_v_t::type::data || inf_rad.count(x)) continue;
            const fp_t& d = m.at(v).at(x);
            if (d < 0) {
                rad_map.erase(x);
                inf_rad.insert(x);
                continue;
            }
            rad_map[x] += d;
        }
    }
    // Select vertex with min radius (this is the center).
    sptr<fpn_v_t> cen = nullptr;
    fp_t min_rad = std::numeric_limits<fp_t>::max();
    for (const auto& [x,r] : rad_map) {
        if (r < min_rad) {
            cen = x;
            min_rad = r;
        }
    }
    return cen;
}

}   // placc
