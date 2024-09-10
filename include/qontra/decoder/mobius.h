/*
 *  author: Suhas Vittal
 *  date:   30 March 2024
 *
 *  NOT CHROMOMBIUS -- This is my own Mobius decoder.
 * */

#ifndef DECODER_MOBIUS_h
#define DECODER_MOBIUS_h

#include "qontra/decoder/restriction.h"

namespace qontra {

namespace gd = graph::decoding;

class MobiusDecoder : public RestrictionDecoder {
public:
    MobiusDecoder(const DetailedStimCircuit&);
    
    // Note that MobiusDecoder's load syndrome ignores the last three operands.
    void load_syndrome(stim::simd_bits_range_ref<SIMD_WIDTH>, int, int, bool) override;
protected:
    // syndrome input is unused.
    std::vector<assign_t> compute_matchings(stim::simd_bits_range_ref<SIMD_WIDTH>) override;

    std::vector<assign_t> compute_matching_on_unified_lattice(void);
    void read_ufl_error_chain(std::vector<assign_t>&, sptr<gd::vertex_t>, sptr<gd::vertex_t>);

    uptr<graph::DecodingGraph> ufl_graph;
    std::unordered_map<sptr<gd::vertex_t>, sptr<gd::vertex_t>> ufl_to_orig_map;
    std::vector<uint64_t> ufl_detectors;
};

}

#endif  // DECODER_MOBIUS_h
