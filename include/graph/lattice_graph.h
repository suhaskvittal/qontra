/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 * */

#ifndef LATTICE_GRAPH_h
#define LATTICE_GRAPH_h

#include "graph/algorithms/distance.h"
#include "graph/graph.h"
#include "graph/io.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace qontra {
namespace graph {

namespace lattice {

struct vertex_t : base::vertex_t {
    uint8_t qubit_type;

    const static uint8_t DATA =     0x0;
    const static uint8_t XPARITY =  0x2;
    const static uint8_t ZPARITY =  0x3;
};

struct edge_t : base::edge_t {
};

}   // lattice

#define __LatticeGraphParent    Graph<lattice::vertex_t, lattice::edge_t>

class LatticeGraph : public __LatticeGraphParent {
public:
    LatticeGraph(void)
        :Graph(),
        meas_order()
    {}

    LatticeGraph(const LatticeGraph& other)
        :Graph(other),
        meas_order(other.meas_order)
    {}

    typedef struct {    // Each pair of vertices has this entry where
                        //  (1) error_chain corresponds to the data qubits
                        //      on the path between the two qubits. This is
                        //      more useful for pairs of parity qubits.
                        //  (2) physical_path is the same as above but also
                        //      includes parity qubits.
        std::vector<lattice::vertex_t*> error_chain;
        std::vector<lattice::vertex_t*> physical_path;
    } matrix_entry_t;

    matrix_entry_t
    get_path_data(lattice::vertex_t* v1, lattice::vertex_t* v2) {
        update_state();
        return distance_matrix[v1][v2];
    }

    lattice::vertex_t*  get_qubit_at_meas_t(uint t) { return meas_order[t]; }

    void    add_to_meas_order(lattice::vertex_t* v) { meas_order.push_back(v); }
protected:
    bool    update_state(void) override;
private:
    void    build_distance_matrix(void);

    std::vector<lattice::vertex_t*> meas_order;
    distance::DistanceMatrix<lattice::vertex_t, matrix_entry_t> distance_matrix;
};

namespace io {

// Function for io callback.
//
// The lattice graph file should be formatted as follows:
//  Data Decl:      D,<data-qubit-1>,<data-qubit-2>,...
//      i.e.            D,0,1,2,3,4
//  Parity Decl:    P<Z,X>,<parity-qubit-1>,<parity-qubit-2>
//      i.e.            PZ,17,18,19,20
//                      PX,21,22,23,24
//                      -- The parity qubits should be listed in the same
//                          order of their measurements.
//  Edge Decl:      E,<edge1-src>,<edge1-dst>,<edge2-src>,<edge2-dst>...
//      i.e.            E,0,17,1,17
//
// The declarations can be split up and can be listed in any order.

void    update_lattice_graph(LatticeGraph&, std::string);

}   // io

}   // graph
}   // qontra

#endif  // LATTICE_GRAPH_h
