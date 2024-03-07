/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#ifndef QONTRA_DECODING_GRAPH_h
#define QONTRA_DECODING_GRAPH_h

#include "qontra/graph/decoding_graph/edge_class.h"
#include "qontra/graph/decoding_graph/structures.h"

#include "qontra/ext/stim.h"
#include "qontra/graph/algorithms/distance.h"

#include <vtils/two_level_map.h>

#include <utility>

namespace qontra {
namespace graph {

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

    sptr<decoding::vertex_t> make_vertex(uint64_t) const override;

    sptr<decoding::hyperedge_t> get_best_shared_edge(std::vector<sptr<decoding::vertex_t>>);

    void immediately_initialize_distances_for(int, int);

    error_chain_t get_error_chain(
            uint64_t,
            uint64_t,
            int=COLOR_ANY,
            int=COLOR_ANY,
            bool force_unflagged=false);

    error_chain_t get_error_chain(
            sptr<decoding::vertex_t>,
            sptr<decoding::vertex_t>,
            int=COLOR_ANY,
            int=COLOR_ANY,
            bool force_unflagged=false);

    std::vector<sptr<decoding::vertex_t>>
        get_complementary_boundaries_to(std::vector<sptr<decoding::vertex_t>>);
    
    void    activate_flags(const std::vector<uint64_t>& all_detectors);
    void    deactivate_flags(void);

    sptr<decoding::vertex_t> get_boundary_vertex(int color);
    sptr<decoding::hyperedge_t> get_base_edge(sptr<decoding::hyperedge_t>);

    std::vector<sptr<decoding::hyperedge_t>> get_flag_edges(void);
    std::vector<sptr<decoding::hyperedge_t>> get_all_edges(void);

    sptr<decoding::hyperedge_t> get_best_edge_from_class_of(sptr<decoding::hyperedge_t>);
    sptr<decoding::hyperedge_t> get_best_nod_edge(bool require_exact_match=false);

    std::map<sptr<decoding::hyperedge_t>, sptr<decoding::hyperedge_t>> get_best_rep_map(void);
    EdgeClass get_edge_class(sptr<decoding::hyperedge_t>);

    poly_t  get_error_polynomial(void);
    fp_t    get_expected_errors(void);

    const int number_of_colors;
protected:
    bool    update_state(void) override;
private:
    sptr<decoding::vertex_t>    make_and_add_vertex_(uint64_t, const DetailedStimCircuit&);

    void resolve_edges(const std::vector<sptr<decoding::hyperedge_t>>&, size_t flips_per_error);

    sptr<decoding::hyperedge_t> get_best_flag_edge(std::vector<sptr<decoding::hyperedge_t>>);

    // If to is not null, then Dijkstra's will terminate upon finding to.
    void dijkstra_(int, int, sptr<decoding::vertex_t> from, sptr<decoding::vertex_t> to=nullptr);
    void make_dijkstra_graph(int, int);
    void build_error_polynomial(void);

    void update_paths(
            uptr<DijkstraGraph>&,
            DistanceMatrix<decoding::vertex_t, error_chain_t>&,
            sptr<decoding::vertex_t> from,
            std::vector<sptr<decoding::vertex_t>> to_list,
            const std::map<sptr<decoding::vertex_t>, fp_t>& dist,
            const std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>& pred);

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

    vtils::TwoLevelMap<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>, fp_t>
        base_probability_map;

    std::set<uint64_t> active_flags;
    std::set<uint64_t> flag_detectors;
    std::vector<EdgeClass> edge_classes;

    // Maps flags to classes that have associated flag edges.
    std::map<uint64_t, std::vector<EdgeClass>> flag_class_map;
    // Maps edges to their containing class.
    std::map<sptr<decoding::hyperedge_t>, EdgeClass> edge_class_map;
    // nod_edges = no detector edges. These are unique flag edges that should be used when
    // flags are active, but no detectors are observed.
    std::vector<sptr<decoding::hyperedge_t>> nod_edges;
    std::vector<sptr<decoding::hyperedge_t>> all_edges;
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
