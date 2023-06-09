/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include "defs.h"
#include "graph/graph.h"

#include <stim.h>

#include <algorithm>
#include <array>

#include <math.h>

namespace qontra {
namespace graph {

const uint BOUNDARY_INDEX = std::numeric_limits<uint>::max();

namespace decoding {

struct vertex_t : base::vertex_t {
    std::array<fp_t, 3> coords;
};

struct edge_t : base::edge_t {
    fp_t            edge_weight;
    fp_t            error_probability;
    std::set<uint>  frames;
};

}   // decoding

#define __DecodingGraphParent   Graph<decoding::vertex_t, decoding::edge_t>

class DecodingGraph : public __DecodingGraphParent {
public:
    DecodingGraph()
        :Graph(), distance_matrix(), error_polynomial(), expected_errors()
    {
        std::array<fp_t, 3> boundary_coords;
        boundary_coords.fill(-1);

        decoding::vertex_t* boundary = new decoding::vertex_t;
        boundary->id = BOUNDARY_INDEX;
        boundary->coords = boundary_coords;
        add_vertex(boundary);
    }

    DecodingGraph(const DecodingGraph& other)
        :Graph(other),
        distance_matrix(other.distance_matrix),
        error_polynomial(other.error_polynomial),
        expected_errors(other.expected_errors)
    {}

    typedef struct {    // Each pair of vertices has this entry, and each entry
                        // corresponds to an error chain.
        uint32_t        chain_length;   // Length of error chain (or shortest path)
        fp_t            probability;    // Probability of error chain
        fp_t            weight;         // Weight used for MWPM decoding
        std::set<uint>  frame_changes;  // Pauli frames that are changed by the chain
    } matrix_entry_t;

    matrix_entry_t
    get_error_chain_data(decoding::vertex_t* v1, decoding::vertex_t* v2) {
        update_state(); 
        return distance_matrix[v1][v2];
    }

    // We can represent the number of errors as a polynomial where the coefficient
    // of xk corresponds to the probability of having k errors.
    //
    // We can compute this polynomial using generating functions.
    //
    // The below functions return that polynomial and the expected number of errors.

    poly_t  get_error_polynomial(void) { update_state(); return error_polynomial; }
    fp_t    get_expected_errors(void) { update_state(); return expected_errors; }
protected:
    bool    update_state(void) override;
private:
    void    build_distance_matrix(void);
    void    build_error_polynomial(void);

    dijkstra::DistanceMatrix<decoding::vertex_t, matrix_entry_t>    distance_matrix;

    poly_t  error_polynomial;
    fp_t    expected_errors;
};

DecodingGraph
to_decoding_graph(const stim::Circuit&);

}   // graph
}   // qontra

#endif
