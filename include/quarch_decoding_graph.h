/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QUARCH_DECODING_GRAPH_h
#define QUARCH_DECODING_GRAPH_h

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <tuple>
#include <map>
#include <functional>

#include <math.h>

#include <stim.h>

#include "quarch_defs.h"

#define N_COORD 100

struct DecodingVertex {
    std::array<fp_t, N_COORD> coord;
    uint detector;
};

struct DecodingEdge {
    std::pair<uint, uint> detectors;
    fp_t edge_weight;
    fp_t error_probability;
    std::set<uint> frames;
};

typedef boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            DecodingVertex,
            DecodingEdge> decoding_graph_base;

#define BOUNDARY_INDEX ((uint)-1)

class DecodingGraph {
public:
    DecodingGraph();
    DecodingGraph(const DecodingGraph&);
    DecodingGraph(DecodingGraph&&);

    void add_detector(uint, std::array<fp_t, N_COORD>& coord);
    void add_edge(uint det1, uint det2, fp_t weight,
            fp_t e_prob, std::set<uint>& frames);
    uint get(uint det_id);

    decoding_graph_base base;
private:
    std::map<uint, uint> detector_to_index;
    std::array<fp_t, N_COORD> boundary_coord;
};

DecodingGraph
to_decoding_graph(const stim::Circuit&);

typedef std::function<void(fp_t, std::vector<uint>, std::set<uint>)>
    error_callback_f;
typedef std::function<void(uint, std::array<fp_t, N_COORD>)>
    detector_callback_f;

void
_read_detector_error_model(const stim::DetectorErrorModel&, 
        uint n_iter, uint& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_f, detector_callback_f);

#endif
