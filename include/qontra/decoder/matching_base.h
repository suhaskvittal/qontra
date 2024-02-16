/*
 *  author: Suhas Vittal
 *  date:   16 February 2024
 * */

#ifndef QONTRA_DECODER_MATCHING_BASE_h
#define QONTRA_DECODER_MATCHING_BASE_h

#include "qontra/decoder.h"
#include "qontra/graph/decoding_graph.h"

namespace qontra {

class MatchingBase : public Decoder {
public:
    MatchingBase(const DetailedStimCircuit& circuit, int flips_per_error)
        :Decoder(circuit),
        decoding_graph(std::make_unique<graph::DecodingGraph>(circuit, flips_per_error)),
        detectors(),
        flags()
    {}

    void load_syndrome(stim::simd_bits_range_ref<SIMD_WIDTH>, int=graph::COLOR_ANY, int=graph::COLOR_ANY);

    std::vector<Decoder::assign_t>
        compute_matching(int=graph::COLOR_ANY, int=graph::COLOR_ANY, bool split_thru_boundary_match=false);
protected:
    uptr<graph::DecodingGraph> decoding_graph;

    std::vector<uint64_t>   detectors;
    std::vector<uint64_t>   flags;
};

}   // qontra
 
#endif  // QONTRA_DECODER_MATCHING_BASE_h
