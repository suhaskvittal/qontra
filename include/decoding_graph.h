/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include "defs.h"

#include <stim.h>

#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <math.h>

#define N_COORD 100

namespace qrc {

struct DijkstraResult {
    std::vector<uint> path;
    fp_t distance;
};

typedef std::map<std::pair<uint, uint>, DijkstraResult> PathTable;

#define BOUNDARY_INDEX ((uint)-1)

class DecodingGraph {
public:
    DecodingGraph();

    struct Vertex {
        int32_t id;
        std::array<fp_t, N_COORD> coord;
        uint detector;

        bool operator==(const Vertex& other) const {
            return id == other.id; 
        }

        bool operator <(const Vertex& other) const {
            return id < other.id; 
        }

        bool operator!=(const Vertex& other) const {
            return !(*this == other); 
        }
    };

    struct Edge {
        int32_t id;
        std::pair<uint, uint> detectors;
        fp_t edge_weight;
        fp_t error_probability;
        std::set<uint> frames;

        bool operator==(const Edge& other) const {
            return id == other.id; 
        }

        bool operator <(const Edge& other) const {
            return id < other.id; 
        }

        bool operator!=(const Edge& other) const {
            return !(*this == other); 
        }
    };

    uint count_detectors(void);

    void add_detector(uint, std::array<fp_t, N_COORD>& coord);
    void add_edge(uint det1, uint det2, fp_t weight,
            fp_t e_prob, std::set<uint>& frames);

    void remove_vertex(const Vertex&);
    void remove_edge(const Edge&);

    Vertex get_vertex(uint det_id);
    Edge get_edge(uint, uint);
   
    std::vector<Vertex> vertices(void);
    std::vector<Vertex> adjacency_list(const Vertex&);
private:
    std::map<uint, Vertex> detector_to_vertex;
    std::map<std::pair<Vertex, Vertex>, Edge> vertices_to_edge;
    std::array<fp_t, N_COORD> boundary_coord;

    std::vector<Vertex> vertex_list;
    std::map<Vertex, std::vector<Vertex>> adjacency_matrix;
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
_dijkstra(DecodingGraph&, const DecodingGraph::Vertex&, 
        std::map<DecodingGraph::Vertex, fp_t>& distances,
        std::map<DecodingGraph::Vertex, DecodingGraph::Vertex>& predecessors);
void
_update_path_table(PathTable&,
        DecodingGraph&, 
        const DecodingGraph::Vertex&,
        const DecodingGraph::Vertex&, 
        std::map<DecodingGraph::Vertex, fp_t>& distances,
        std::map<DecodingGraph::Vertex, DecodingGraph::Vertex>& predecessors);

}  // qrc

#endif
