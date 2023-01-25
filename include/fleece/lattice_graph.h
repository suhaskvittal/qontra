/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_LATTICE_GRAPH_h
#define FLEECE_LATTICE_GRAPH_h

#include <stim.h>

#include "decoding_graph.h"
#include "defs.h"
#include "graph/dijkstra.h"

#include <deque>
#include <map>
#include <string>
#include <vector>

namespace qrc {
namespace fleece {

#define LATTICE_BOUNDARY -1
#define LATTICE_BOUNDARY_DETECTOR -1

class LatticeGraph {
public:
    LatticeGraph();
    ~LatticeGraph();

    struct Vertex {
        int32_t qubit;
        bool is_data;
        std::vector<int32_t> detectors;
        std::vector<uint32_t> measurement_times;

        Vertex()
            :qubit(-1),
            is_data(false),
            detectors(),
            measurement_times()
        {}

        Vertex(int32_t q, bool is_data)
            :qubit(q),
            is_data(is_data),
            detectors(),
            measurement_times()
        {}

        Vertex(const Vertex& other)
            :qubit(other.qubit), 
            is_data(other.is_data), 
            detectors(other.detectors),
            measurement_times(other.measurement_times)
        {}

        bool operator==(const Vertex& other) const {
            return qubit == other.qubit;
        }

        bool operator<(const Vertex& other) const {
            return qubit < other.qubit;
        }
    };

    // Adding returns false if the qubit already exists
    // Data may be updated though.
    bool add_qubit(int32_t, uint8_t is_data, int32_t detector=-1, 
            int32_t meas_time=-1, bool force_record=false);
    bool add_coupling(int32_t, int32_t);
    bool add_coupling(Vertex*, Vertex*);

    Vertex* get_vertex_by_qubit(int32_t);
    Vertex* get_vertex_by_detector(int32_t);

    std::vector<Vertex*> vertices(void);

    std::vector<Vertex*> adjacency_list(int32_t);
    std::vector<Vertex*> adjacency_list(Vertex*);

    std::vector<Vertex*> get_common_neighbors(Vertex*, Vertex*);
private:
    std::vector<Vertex*> vertex_list;
    std::map<int32_t, Vertex*> qubit_to_vertex;
    std::map<int32_t, Vertex*> detector_to_vertex;
    std::map<Vertex*, std::vector<Vertex*>> adjacency_matrix;

    friend LatticeGraph to_lattice_graph(const stim::Circuit&);
};

LatticeGraph 
to_lattice_graph(const stim::Circuit&);
graph::PathTable<LatticeGraph::Vertex>
compute_path_table(LatticeGraph&);

}  // fleece
}  // qrc


#endif
