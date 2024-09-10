/*
 *  author: Suhas Vittal
 *  date:   23 March 2024
 * */

#ifndef PROTEAN_YIELD_SIM_h
#define PROTEAN_YIELD_SIM_h

#include "qontra/protean/network.h"

#include <unordered_map>

namespace qontra {
namespace protean {

const fp_t ANH = -340e6;

class YieldSimulator {
public:
    YieldSimulator(PhysicalNetwork*);

    fp_t est_mean_collisions(fp_t fab_precision,
                                uint64_t trials,
                                uint64_t seed=0,
                                sptr<net::phys_vertex_t> =nullptr);
    // freq_list is a list of candidate frequencies. It is assumed to be sorted.
    void assign(fp_t fab_precision, const std::vector<fp_t>& freq_list);

    std::unordered_map<sptr<net::phys_vertex_t>, fp_t> freq_map;
private:
    size_t count_violations(
            sptr<net::phys_vertex_t>, const std::unordered_map<sptr<net::phys_vertex_t>, fp_t>&);

    size_t count_collisions(const std::unordered_map<sptr<net::phys_vertex_t>, fp_t>&, 
                            const std::vector<sptr<net::phys_vertex_t>>& vertices);
    size_t compute_center_score(sptr<net::phys_vertex_t>);

    fp_t local_sim(sptr<net::phys_vertex_t>, fp_t fab_precision, uint64_t trials);

    PhysicalNetwork* network;
    sptr<net::phys_vertex_t> vcenter;
};

}   // protean
}   // qontra

#endif  // PROTEAN_YIELD_SIM_h
