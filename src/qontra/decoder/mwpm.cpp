/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "qontra/decoder/mwpm.h"

#include <PerfectMatching.h>

namespace qontra {

using namespace graph;
using namespace decoding;

MWPMDecoder::MWPMDecoder(const DetailedStimCircuit& circuit)
    :MatchingBase(circuit, 2)
{}

Decoder::result_t
MWPMDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();
    stim::simd_bits<SIMD_WIDTH> corr(n_obs);

    load_syndrome(syndrome);

    auto _detectors = detectors;
    auto _flags = flags;

    if (detectors.empty()) return ret_no_detectors();

    timer.clk_start();
    std::vector<assign_t> assignments = compute_matching();
    fp_t t = static_cast<fp_t>(timer.clk_end());
    
#ifdef DECODER_PERF
    fp_t _t;
    timer.clk_start();
#endif

    for (const assign_t& m : assignments) {
        for (size_t i = 1; i < m.path.size(); i++) {
            sptr<vertex_t> v = m.path[i-1],
                            w = m.path[i];
            sptr<hyperedge_t> e = decoding_graph->get_edge({v, w}); 
            if (e == nullptr) {
                // First try and get the flag edge.
                e = get_flag_edge_for({v, w});
            }
            // If it is still not found, then just declare an error and move on.
            if (e == nullptr) {
                continue;
            }
            e = decoding_graph->get_best_edge_from_class_of(e);
            for (uint64_t f : e->frames) corr[f] ^= 1;
        }
    }

#ifdef DECODER_PERF
    _t = timer.clk_end();
    std::cout << "[ MWPMDecoder ] took " << _t*1e-9 << "s to retrieve correction from matching" << std::endl;
#endif

    return { t, corr };
}

}   // qontra
