/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QONTRA_DECODING_GRAPH_h
#define QONTRA_DECODING_GRAPH_h

#include "qontra/defs.h"
#include "qontra/graph/algorithms/distance.h"
#include "qontra/graph.h"

#include <stim.h>

namespace qontra {
namespace graph {

const uint64_t BOUNDARY_INDEX = std::numeric_limits<uint32_t>::max();

const uint64_t RED_BOUNDARY_INDEX = (1L << 48) | BOUNDARY_INDEX;
const uint64_t BLUE_BOUNDARY_INDEX = (1L << 49) | BOUNDARY_INDEX;
const uint64_t GREEN_BOUNDARY_INDEX = (1L << 50) | BOUNDARY_INDEX;

namespace decoding {

#define N_COORD    16

struct vertex_t : base::vertex_t {
    std::array<fp_t, N_COORD>    coords;
};

struct edge_t : base::edge_t {
    fp_t                edge_weight;
    fp_t                error_probability;
    std::set<uint64_t>  frames;
};

struct colored_vertex_t : vertex_t {
    std::string color;
};

struct colored_edge_t : edge_t {
};

}   // decoding

bool    is_boundary(sptr<decoding::vertex_t>);
bool    is_colored_boundary(sptr<decoding::vertex_t>);

std::string int_to_color(int x);
int         color_to_int(std::string x);

#define __DecodingGraphParent   Graph<decoding::vertex_t, decoding::edge_t>

class DecodingGraph : public __DecodingGraphParent {
public:
    enum class Mode {
        NORMAL,     // No adjustments. Default option.
        LOW_MEMORY, // Here, we will not generate the distance matrix.
        DO_NOT_BUILD,   // Hidden option for decoders to avoid building the decoding graph.
    };

    DecodingGraph(Mode=Mode::NORMAL);
    DecodingGraph(const DecodingGraph& other);

    DecodingGraph& operator=(const DecodingGraph& other);

    typedef struct {    // Each pair of vertices has this entry, and each entry
                        // corresponds to an error chain.
        size_t  chain_length;   // Length of error chain (or shortest path)
        fp_t    probability;    // Probability of error chain
        fp_t    weight;         // Weight used for MWPM decoding

        std::set<uint64_t>  frame_changes;  // Pauli frames that are changed by the chain

        std::vector<sptr<decoding::vertex_t>> error_chain;
        bool error_chain_runs_through_boundary;
    } matrix_entry_t;

    matrix_entry_t get_error_chain_data(sptr<decoding::vertex_t>, sptr<decoding::vertex_t>);

    // We can represent the number of errors as a polynomial where the coefficient
    // of xk corresponds to the probability of having k errors.
    //
    // We can compute this polynomial using generating functions.
    //
    // The below functions return that polynomial and the expected number of errors.

    poly_t  get_error_polynomial(void);
    fp_t    get_expected_errors(void);
protected:
    bool    update_state(void) override;
private:
    matrix_entry_t dijkstra_cb(
                        sptr<decoding::vertex_t>,
                        sptr<decoding::vertex_t>,
                        const std::map<sptr<decoding::vertex_t>, fp_t>&,
                        const std::map<sptr<decoding::vertex_t>, sptr<decoding::vertex_t>>&);

    void    build_distance_matrix(void);
    void    build_flagged_decoding_graph(void); // This is a partial copy.
    void    build_error_polynomial(void);

    DistanceMatrix<decoding::vertex_t, matrix_entry_t>  distance_matrix;

    poly_t  error_polynomial;
    fp_t    expected_errors;

    Mode mode;
};

// This is standard method of building decoding graphs, where each error is assumed
// to flip exactly two detectors. This works for codes like the surface code.
//
// For codes such as the color code where an error can flip 3 detectors, a custom 
// method must be used. We leave this at the discretion of the user, but we note
// this can be easily performed by leveraging the metadata for each detector
// in the Stim circuit (i.e. instead of coordinates, use colors and have your own
// method use that).

DecodingGraph
to_decoding_graph(const stim::Circuit&, DecodingGraph::Mode=DecodingGraph::Mode::NORMAL);

}   // graph
}   // qontra

#include "decoding_graph.inl"

#endif  // QONTRA_DECODING_GRAPH_h
