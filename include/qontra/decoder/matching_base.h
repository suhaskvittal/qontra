/*
 *  author: Suhas Vittal
 *  date:   16 February 2024
 * */

#ifndef QONTRA_DECODER_MATCHING_BASE_h
#define QONTRA_DECODER_MATCHING_BASE_h

#include "qontra/decoder.h"
#include "qontra/graph/decoding_graph.h"

namespace qontra {

namespace gd = graph::decoding;

struct assign_t {
    sptr<gd::vertex_t> v = nullptr;
    sptr<gd::vertex_t> w = nullptr;
    int c1;
    int c2;
    std::vector<sptr<gd::vertex_t>> path;

    inline bool operator==(const assign_t& other) const {
        return as_tuple() == other.as_tuple();
    }

    inline bool operator<(const assign_t& other) const {
        return as_tuple() < other.as_tuple();
    }
private:
    std::tuple<sptr<gd::vertex_t>, sptr<gd::vertex_t>, int, int, std::vector<sptr<gd::vertex_t>>>
        as_tuple(void) const { return std::make_tuple(v, w, c1, c2, path); }
};

class MatchingBase : public Decoder {
public:

    MatchingBase(const DetailedStimCircuit&, int flips_per_error);

    void load_syndrome(
            stim::simd_bits_range_ref<SIMD_WIDTH>,
            int=graph::COLOR_ANY,
            int=graph::COLOR_ANY,
            bool recompute_flags=true);

    stim::simd_bits<SIMD_WIDTH> get_base_corr(void);
    Decoder::result_t ret_no_detectors(void);

    std::vector<assign_t> compute_matching(
            int=graph::COLOR_ANY, int=graph::COLOR_ANY, bool split_thru_boundary_match=false);

    sptr<graph::decoding::hyperedge_t> get_flag_edge_for(std::vector<sptr<gd::vertex_t>>);
protected:
    graph::error_chain_t expand_error_chain(
            uint64_t, uint64_t, int=graph::COLOR_ANY, int=graph::COLOR_ANY);
    graph::error_chain_t expand_error_chain(
            sptr<gd::vertex_t>, sptr<gd::vertex_t>, int=graph::COLOR_ANY, int=graph::COLOR_ANY);
    void update_assignments_wrt_boundary(std::vector<assign_t>&, uint64_t, uint64_t, int, int);

    uptr<graph::DecodingGraph> decoding_graph;

    std::vector<uint64_t>   detectors;
    std::vector<uint64_t>   flags;
private:
    std::vector<sptr<gd::hyperedge_t>> flag_edges;
};

}   // qontra
 
#endif  // QONTRA_DECODER_MATCHING_BASE_h
