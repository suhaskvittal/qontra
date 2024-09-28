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
    assign_t()
        :v(nullptr),
        w(nullptr),
        c1(graph::COLOR_ANY),
        c2(graph::COLOR_ANY),
        flag_edges(),
        path()
    {}

    ~assign_t() {}

    sptr<gd::vertex_t> v;
    sptr<gd::vertex_t> w;
    int c1;
    int c2;
    std::vector<sptr<gd::hyperedge_t>> flag_edges;
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
    MatchingBase(
        const DetailedStimCircuit&,
        int flips_per_error,
        bool reweigh_for_detectors=false);

    virtual void load_syndrome(
            stim::simd_bits_range_ref<SIMD_WIDTH>, int=graph::COLOR_ANY, int=graph::COLOR_ANY,
            bool recompute_flags=true);

    virtual std::vector<assign_t> compute_matching(int=graph::COLOR_ANY, int=graph::COLOR_ANY);

    Decoder::result_t ret_no_detectors(void);
    sptr<gd::hyperedge_t> get_flag_edge_for(std::vector<sptr<gd::vertex_t>>);
    sptr<gd::hyperedge_t> get_flag_edge_for(std::vector<sptr<gd::vertex_t>>, const std::vector<sptr<gd::hyperedge_t>>&);
protected:
    void identify_flag_edges_in_path(assign_t&);

    uptr<graph::DecodingGraph> decoding_graph;

    std::vector<uint64_t>   detectors;
    std::vector<uint64_t>   flags;
private:
    std::vector<sptr<gd::hyperedge_t>> flag_edges;
};

}   // qontra

#endif  // QONTRA_DECODER_MATCHING_BASE_h
