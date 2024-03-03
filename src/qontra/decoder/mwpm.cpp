/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "qontra/decoder/mwpm.h"

#include <PerfectMatching.h>

namespace qontra {

using namespace graph;
using namespace decoding;

Decoder::result_t
MWPMDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();
    stim::simd_bits<SIMD_WIDTH> corr(n_obs);

    load_syndrome(syndrome);

    auto _detectors = detectors;
    auto _flags = flags;

    /*
    std::cout << "syndrome: D[";
    for (uint64_t d : _detectors) std::cout << " " << d;
    std::cout << " ], F[";
    for (uint64_t f : _flags) std::cout << " " << f;
    std::cout << " ]" << std::endl;
    */

    if (detectors.empty()) return ret_no_detectors();

    corr ^= get_base_corr();

    timer.clk_start();
    std::vector<Decoder::assign_t> assignments = compute_matching();
    fp_t t = static_cast<fp_t>(timer.clk_end());

    for (auto& match : assignments) {
        uint64_t di = std::get<0>(match),
                 dj = std::get<1>(match);
        error_chain_t ec = decoding_graph->get(di, dj);
        for (size_t i = 1; i < ec.path.size(); i++) {
            sptr<vertex_t> vi = ec.path[i-1],
                            vj = ec.path[i];
            sptr<hyperedge_t> e = decoding_graph->get_edge({vi, vj}); 
            if (e == nullptr) {
                // First try and get the flag edge.
                e = get_flag_edge_for({vi, vj});
            }
            // If it is still not found, then just declare an error and move on.
            if (e == nullptr) {
                continue;
            }
            e = decoding_graph->get_best_edge_from_class_of(e);
            for (uint64_t f : e->frames) corr[f] ^= 1;
        }
    }

    return { t, corr, assignments };
}

}   // qontra
