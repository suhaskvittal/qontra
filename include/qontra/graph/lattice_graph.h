/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 * */

#ifndef QONTRA_LATTICE_GRAPH_h
#define QONTRA_LATTICE_GRAPH_h

#include "qontra/graph/algorithms/distance.h"
#include "qontra/graph/graph.h"

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
    LatticeGraph(void);
    LatticeGraph(const LatticeGraph&);

    typedef struct {    // Each pair of vertices has this entry where
                        //  (1) error_chain corresponds to the data qubits
                        //      on the path between the two qubits. This is
                        //      more useful for pairs of parity qubits.
                        //  (2) physical_path is the same as above but also
                        //      includes parity qubits.
        std::vector<sptr<lattice::vertex_t>> error_chain;
        std::vector<sptr<lattice::vertex_t>> physical_path;
    } matrix_entry_t;

    matrix_entry_t get_path_data(sptr<lattice::vertex_t>, sptr<lattice::vertex_t>);

    std::vector<std::vector<sptr<lattice::vertex_t>>> x_obs_list;
    std::vector<std::vector<sptr<lattice::vertex_t>>> z_obs_list;
protected:
    bool    update_state(void) override;
private:
    void    build_distance_matrix(void);

    DistanceMatrix<lattice::vertex_t, matrix_entry_t> distance_matrix;
};

}   // graph
}   // qontra

#include "lattice_graph.inl"

#endif  // LATTICE_GRAPH_h
