/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include "defs.h"
#include "graph/graph.h"
#include "graph/algorithms/distance.h"

#include <stim.h>

#include <algorithm>
#include <array>

#include <math.h>
namespace qontra {
namespace graph {

const uint BOUNDARY_INDEX = std::numeric_limits<uint>::max();

const uint64_t RED_BOUNDARY_INDEX = (1L << 48) | BOUNDARY_INDEX;
const uint64_t BLUE_BOUNDARY_INDEX = (1L << 49) | BOUNDARY_INDEX;
const uint64_t GREEN_BOUNDARY_INDEX = (1L << 50) | BOUNDARY_INDEX;

namespace decoding {

#define N_COORD    16

struct vertex_t : base::vertex_t {
    std::array<fp_t, N_COORD>    coords;
};

struct edge_t : base::edge_t {
    fp_t            edge_weight;
    fp_t            error_probability;
    std::set<uint>  frames;
};

struct colored_vertex_t : vertex_t {
    std::string color;
};

struct colored_edge_t : edge_t {
};

}   // decoding

#define __DecodingGraphParent   Graph<decoding::vertex_t, decoding::edge_t>

class DecodingGraph : public __DecodingGraphParent {
public:
    enum class Mode {
        NORMAL,     // No adjustments. Default option.
        LOW_MEMORY, // Here, we will not generate the distance matrix.
        DO_NOT_BUILD,   // Hidden option for decoders to avoid building the decoding graph.
    };

    DecodingGraph(Mode m=Mode::NORMAL)
        :Graph(), distance_matrix(), error_polynomial(), expected_errors(), mode(m)
    {
        std::array<fp_t, N_COORD> boundary_coords;
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
        expected_errors(other.expected_errors),
        mode(other.mode)
    {}

    typedef struct {    // Each pair of vertices has this entry, and each entry
                        // corresponds to an error chain.
        uint32_t        chain_length;   // Length of error chain (or shortest path)
        fp_t            probability;    // Probability of error chain
        fp_t            weight;         // Weight used for MWPM decoding
        std::set<uint>  frame_changes;  // Pauli frames that are changed by the chain

        std::vector<decoding::vertex_t*>    error_chain;
        bool                                error_chain_runs_through_boundary;
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

    distance::DistanceMatrix<decoding::vertex_t, matrix_entry_t>    distance_matrix;

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

// Colored decoding graphs are a collection of decoding graphs, such that each decoding
// graph restricts the decoding graph to two colors. Now, each detector in the graph corresponds
// to detectors in multiple decoding graphs, and we must track such assignments alongside
// other aspects such as adjacency.

#define __ColoredDecodingGraphParent    Graph<decoding::colored_vertex_t, decoding::colored_edge_t>

typedef std::tuple<decoding::colored_vertex_t*,
                    decoding::colored_vertex_t*,
                    decoding::colored_vertex_t*> face_t;
face_t
make_face(decoding::colored_vertex_t*, decoding::colored_vertex_t*, decoding::colored_vertex_t*);

inline std::string int_to_color(int x) {
    return std::to_string("rgb"[x]);
}

class ColoredDecodingGraph : public __ColoredDecodingGraphParent {
public:
    ColoredDecodingGraph(DecodingGraph::Mode mode=DecodingGraph::Mode::NORMAL);

    bool    add_vertex(decoding::colored_vertex_t*) override;
    bool    add_edge(decoding::colored_edge_t*) override;

    void    delete_vertex(decoding::colored_vertex_t*) override;
    void    delete_edge(decoding::colored_edge_t*) override;

    std::set<face_t>    get_all_incident_faces(decoding::colored_vertex_t*);
    stim::simd_bits     get_correction_for_face(face_t);

    DecodingGraph& operator[](const std::string& cc) {
        return restricted_graphs.at(restricted_color_map.at(cc));
    }

    DecodingGraph& operator[](const char* cc) {
        return (*this)[std::string(cc)];
    }
private:
    std::map<std::string, int>      restricted_color_map;
    std::array<DecodingGraph, 3>    restricted_graphs;
};

ColoredDecodingGraph
to_colored_decoding_graph(const stim::Circuit&, DecodingGraph::Mode=DecodingGraph::Mode::NORMAL);

}   // graph
}   // qontra

#endif
