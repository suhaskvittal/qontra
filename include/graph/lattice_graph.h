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
    enum class type { data, xparity, zparity, flag };

    type qubit_type;
};

struct edge_t : base::edge_t {
    uint cx_time;
};

}   // lattice

#define __LatticeGraphParent    Graph<lattice::vertex_t, lattice::edge_t>

class LatticeGraph : public __LatticeGraphParent {
public:
    LatticeGraph(void)
        :Graph()
    {}

    LatticeGraph(const LatticeGraph& other)
        :Graph(other),
        x_obs_list(other.x_obs_list),
        z_obs_list(other.z_obs_list)
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

    std::vector<std::vector<lattice::vertex_t*>> z_obs_list;
    std::vector<std::vector<lattice::vertex_t*>> x_obs_list;
protected:
    bool    update_state(void) override;
private:
    void    build_distance_matrix(void);

    distance::DistanceMatrix<lattice::vertex_t, matrix_entry_t> distance_matrix;
};

}   // graph
}   // qontra

#endif  // LATTICE_GRAPH_h
