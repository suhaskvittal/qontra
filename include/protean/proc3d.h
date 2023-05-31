/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#ifndef PROTEAN_PROC3D_h
#define PROTEAN_PROC3D_h

#include "defs.h"
#include "graph/coupling_graph.h"

#include <limits>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

namespace proc3d {

struct vertex_t : coupling::vertex_t {
    uint processor_layer;

    bool is_tsv_junction(void) const { return processor_layer > 0; }
};

struct edge_t : coupling::edge_t {
    bool is_vertical;
};

}   // proc3d

class Processor3D : public Graph<proc3d::vertex_t, proc3d::edge_t> {
public:
    Processor3D(uint max_thickness)   
        :CouplingGraph(false),
        main_processor(true),
        processor_layers(),
        max_thickness(max_thickness),
        vertex_to_tsv_junctions()
    {}
    Processor3D(const Processor3D& other)
        :CouplingGraph(other),
        main_processor(other.main_processor),
        processor_layers(other.processor_layers),
        max_thickness(other.max_thickness),
        vertex_to_tsv_junctions(other.vertex_to_tsv_junctions)
    {}

    bool        add_edge(proc3d::edge_t*, bool is_undirected=true) override;

    // The main processor contains the physical qubits and some couplings. The other
    // processor layers contains couplings enabled by TSVs. Each of the processor layers 
    // is guaranteed to be planar. Note that we only need to track planarity for the main
    // processor because the secondary layers only have vertices that are TSVs, and each
    // of these can only have one connection.
    CouplingMap get_main_processor(void) { return main_processor; }
    CouplingMap get_processor_layer(uint k) { return processor_layer[k]; }
    uint        get_thickness(void) { return processor_layers.size(); }
private:
    CouplingGraph               main_processor;
    std::vector<CouplingGraph>  processor_layers;
    uint                        max_thickness;

    std::map<proc3d::vertex_t*, std::vector<proc3d::vertex_t*>> vertex_to_tsv_junctions;   
                                                                        // Maps a vertex on
                                                                        // the main processor
                                                                        // to the TSV vertices
                                                                        // above it on the
                                                                        // processor layers.
};

}   // protean
}   // qontra

#endif  // PROTEAN_PROC3D_h
