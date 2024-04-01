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
    
    bool runs_through_boundary = false;
    std::vector<sptr<decoding::vertex_t>>   path;
    std::vector<sptr<decoding::vertex_t>>   boundary_vertices;
};

typedef std::vector<fp_t> poly_t;
typedef Graph<decoding::vertex_t, decoding::edge_t> DijkstraGraph;

class DecodingGraph : public HyperGraph<decoding::vertex_t, decoding::hyperedge_t> {
public:
    DecodingGraph(const DetailedStimCircuit&, size_t flips_per_error);
    DecodingGraph(DecodingGraph&&) = default;

    DecodingGraph make_unified_lattice(std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>& ufl_map);

    // Makes vertex and also sets the base of the vertex to itself.
    sptr<decoding::vertex_t> make_vertex(uint64_t) const override;
    // Gets shared edge with most similarity to the passed in inputs.
    sptr<decoding::hyperedge_t> get_best_shared_edge(std::vector<sptr<decoding::vertex_t>>);
    // Initializes the distance matrix for the two colors using Floyd-Warshall.
    void init_distances_for(int=COLOR_ANY, int=COLOR_ANY);

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
    
    // activate_detectors is useful for flagged decoding, as doing so will (1) activate
    // any flags in the syndrome, and (2) restrict thee size of the decoding graph
    // when computing Dijkstra's.
    void    activate_detectors(const std::vector<uint64_t>& all_detectors);
    void    activate_detectors(const std::vector<uint64_t>& nonflags, const std::vector<uint64_t>& flags);
    void    deactivate_detectors(void);

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
    DecodingGraph(void);

    bool    update_state(void) override;
private:
    // Add vertices while also adding their bases. This function recursively calls itself until
    // we find a detector whose base is itself.
    sptr<decoding::vertex_t> make_and_add_vertex_(uint64_t, const DetailedStimCircuit&);
    // This function creates equivalence classes and initializes the DecodingGraph accordingly.
    void resolve_edges(const std::vector<sptr<decoding::hyperedge_t>>&, size_t flips_per_error);

    fp_t compute_renorm_factor(std::set<uint64_t> flags={});
    sptr<decoding::hyperedge_t> get_best_flag_edge(std::vector<sptr<decoding::hyperedge_t>>);

    // If to is not null, then Dijkstra's will terminate upon finding to.
    void dijkstra_(int, int, sptr<decoding::vertex_t> from, sptr<decoding::vertex_t> to=nullptr);
    void make_dijkstra_graph(int, int);
    // General path update function used by dijkstra_ and immediately_initalize_distances_for
    void update_paths(
            uptr<DijkstraGraph>&,
            DistanceMatrix<decoding::vertex_t, error_chain_t>&,
            sptr<decoding::vertex_t> from,
            std::vector<sptr<decoding::vertex_t>> to_list,
            const std::map<sptr<decoding::vertex_t>, fp_t>& dist,
            const std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>& pred);

    void build_error_polynomial(void);

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

    std::set<uint64_t> active_detectors;
    std::set<uint64_t> active_flags;
    // flag_detectors is a list of all flag detectors (not just those in a syndrome). This is
    // used to construct active_detectors and active_flags.
    std::set<uint64_t> flag_detectors;
    // List of equivalence classes
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

    fp_t renorm_factor;
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
