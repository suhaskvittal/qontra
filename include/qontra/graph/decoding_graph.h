/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#ifndef QONTRA_DECODING_GRAPH_h
#define QONTRA_DECODING_GRAPH_h

#include "qontra/ext/stim.h"
#include "qontra/hypergraph.h"
#include "qontra/graph/algorithms/distance.h"

#include <limits>
#include <utility>

namespace qontra {
namespace graph {

const int COLOR_ANY = -1;
const uint64_t BOUNDARY_FLAG = std::numeric_limits<uint32_t>::max();

namespace decoding {

struct vertex_t : base::vertex_t {
    int             color = COLOR_ANY;
    sptr<vertex_t>  base = nullptr;
    bool            is_boundary_vertex = false;

    sptr<vertex_t> get_base(void);
};

struct edge_t : base::edge_t {
    fp_t probability;
};

struct hyperedge_t : base::hyperedge_t {
    fp_t                probability;
    std::set<uint64_t>  flags;
    std::set<uint64_t>  frames;
};

}   // decoding

struct error_chain_t {
    size_t  length = 0;
    fp_t    probability = 1.0;
    fp_t    weight = 0.0;
    
    bool runs_through_boundary;
    std::vector<sptr<decoding::vertex_t>>   path;
    std::vector<sptr<decoding::vertex_t>>   boundary_vertices;
};

typedef std::vector<fp_t> poly_t;

typedef Graph<decoding::vertex_t, decoding::edge_t> DijkstraGraph;

class DecodingGraph : public HyperGraph<decoding::vertex_t, decoding::hyperedge_t> {
public:
    DecodingGraph(const DetailedStimCircuit&, size_t flips_per_error);
    DecodingGraph(DecodingGraph&&) = default;

    sptr<decoding::vertex_t>    make_vertex(uint64_t) const override;

    sptr<decoding::vertex_t>    get_boundary_vertex(int color);

    error_chain_t get(uint64_t, uint64_t, bool force_unflagged=false);
    error_chain_t get(sptr<decoding::vertex_t>, sptr<decoding::vertex_t>, bool force_unflagged=false);
    error_chain_t get(int c1, int c2, uint64_t, uint64_t, bool force_unflagged=false);
    error_chain_t get(
            int c1, int c2, sptr<decoding::vertex_t>, sptr<decoding::vertex_t>, bool force_unflagged=false);

    std::vector<sptr<decoding::vertex_t>>
        get_complementary_boundaries_to(std::vector<sptr<decoding::vertex_t>>);
    
    void    activate_flags(const std::vector<uint64_t>& all_detectors);
    void    deactivate_flags(void);

    std::vector<sptr<decoding::hyperedge_t>> get_flag_edges(void);
    std::vector<sptr<decoding::hyperedge_t>> get_flag_singletons(void);

    poly_t  get_error_polynomial(void);
    fp_t    get_expected_errors(void);

    const int number_of_colors;
protected:
    bool    update_state(void) override;
private:
    sptr<decoding::vertex_t>    make_and_add_vertex_(uint64_t, const DetailedStimCircuit&);

    std::vector<sptr<decoding::hyperedge_t>>
        get_flags_from_map(const std::map<uint64_t, std::set<sptr<decoding::hyperedge_t>>>&);

    void    dijkstra_(int, int, sptr<decoding::vertex_t> from);
    void    make_dijkstra_graph(int, int);
    void    build_error_polynomial(void);

    poly_t  error_polynomial;
    fp_t    expected_errors;

    std::map<std::pair<int, int>, uptr<DijkstraGraph>>
        dijkstra_graph_map;
    std::map<std::pair<int, int>, uptr<DijkstraGraph>>
        flagged_dijkstra_graph_map;
    std::map<std::pair<int, int>, DistanceMatrix<decoding::vertex_t, error_chain_t>>
        distance_matrix_map;
    std::map<std::pair<int, int>, DistanceMatrix<decoding::vertex_t, error_chain_t>>
        flagged_distance_matrix_map;

    std::vector<uint64_t> active_flags;

    std::set<uint64_t> flag_detectors;
    std::map<uint64_t, std::set<sptr<decoding::hyperedge_t>>> flag_singleton_map;
    std::map<uint64_t, std::set<sptr<decoding::hyperedge_t>>> flag_edge_map;
    bool flags_are_active;
};

std::vector<int> get_complementary_colors_to(std::vector<int>, int number_of_colors);

uint64_t    get_color_boundary_index(int);
fp_t        compute_weight(fp_t probability);

// DEM_ERROR_FUNC should take in (1) an error probability, (2) a vector of
// uint64_t corresponding to the detectors, and (3) a set of uint64_t
// corresponding to the frame changes.
//
// DEM_DET_FUNC should take in a single parameter: a detector.
template <class DEM_ERROR_FUNC, class DEM_DET_FUNC>
void read_detector_error_model(
        const stim::DetectorErrorModel&,
        size_t n_iter,
        size_t& detector_offset,
        DEM_ERROR_FUNC,
        DEM_DET_FUNC);

}   // graph
}   // qontra

#include "decoding_graph.inl"

#endif  // QONTRA_DECODING_GRAPH_h
