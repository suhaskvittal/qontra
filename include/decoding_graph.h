/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <tuple>
#include <map>
#include <functional>

#include <math.h>

#include <stim.h>

#include "defs.h"

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

struct DijkstraResult {
    std::vector<uint> path;
    fp_t distance;
};

typedef boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            DecodingVertex,
            DecodingEdge> decoding_graph_base;
typedef std::map<std::pair<uint, uint>, DijkstraResult> PathTable;

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

PathTable
compute_path_table(DecodingGraph&);

typedef std::function<void(fp_t, std::vector<uint>, std::set<uint>)>
    error_callback_f;
typedef std::function<void(uint, std::array<fp_t, N_COORD>)>
    detector_callback_f;

void
_read_detector_error_model(const stim::DetectorErrorModel&, 
        uint n_iter, uint& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_f, detector_callback_f);

void
_update_path_table(PathTable&,
        DecodingGraph&, uint, uint, 
        const std::vector<fp_t>& distances,
        const std::vector<uint>& predecessors);

#endif
