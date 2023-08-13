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

#include <lemon/list_graph.h>
#include <lemon/planarity.h>

#include <math.h>

namespace qontra {
namespace protean {

namespace proc3d {

struct vertex_t : graph::coupling::vertex_t {
    uint processor_layer;

    bool is_tsv_junction(void) const { return processor_layer > 0; }
};

struct edge_t : graph::coupling::edge_t {
    bool is_vertical;
    uint processor_layer;
};

}   // proc3d

// The Processor3D's graph structure does not represent the physical
// characteristics of the chip, and instead represents the connectivity
// between different qubits on the chip. However, any methods involving
// edges have been modified to reflect the physical architecture. 
//
// For example, the add_edge method appends an edge in a way that maintains
// the planarity of the main processor (containing the physical qubits) and
// the other layers through which TSVs may pass through. Get edge also does
// so similarly and returns the set of edges which physically connect
// two physical qubits (either one edge (on main processor) or three edges
// (two edges required for TSV)).
#define __Processor3DParent    graph::Graph<proc3d::vertex_t, proc3d::edge_t>

class Processor3D : public __Processor3DParent {
public:
    Processor3D(uint max_thickness=std::numeric_limits<uint>::max())
        :__Processor3DParent(),
        main_processor(true),
        processor_layers(),
        max_thickness(max_thickness), vertex_to_tsv_junctions(),
        coupling_to_edges()
    {
        main_processor.dealloc_on_delete = false;
    }

    Processor3D(const Processor3D& other)
        :__Processor3DParent(other),
        main_processor(other.main_processor),
        processor_layers(other.processor_layers),
        max_thickness(other.max_thickness),
        vertex_to_tsv_junctions(other.vertex_to_tsv_junctions),
        coupling_to_edges(other.coupling_to_edges)
    {}

    ~Processor3D(void) {
        for (auto pair : coupling_to_edges) {
            auto impl = pair.second;
            // Delete physical edges for any virtual edges.
            if (impl.size() > 1) {
                delete impl[0];
                delete impl[2];
            }
        }
    }

    bool    add_vertex(proc3d::vertex_t*) override;
    bool    add_edge(proc3d::edge_t*, bool force_out_of_plane);
    void    delete_vertex(proc3d::vertex_t*) override;
    void    delete_edge(proc3d::edge_t*) override;

    bool    add_edge(proc3d::edge_t* e) override { return add_edge(e, false); }

    void    reallocate_edges(void); // This function greedily tries to add vertical edges
                                    // to the main processor.

    // This method will retrieve the edges corresponding to the physical coupling
    // (i.e. the vertical TSVs are included).
    std::vector<proc3d::edge_t*>    get_physical_edges(proc3d::vertex_t* v1, proc3d::vertex_t* v2)
                { auto e = get_edge(v1, v2); return get_physical_edges(e); }
    std::vector<proc3d::edge_t*>    get_physical_edges(proc3d::edge_t* e)
                { return coupling_to_edges[e]; }

    bool    has_complex_coupling(proc3d::vertex_t* v1, proc3d::vertex_t* v2)
                { auto e = get_edge(v1, v2); return has_complex_coupling(e); }
    bool    has_complex_coupling(proc3d::edge_t* e) { return coupling_to_edges[e].size() > 1; }

    void    force_out_of_plane(proc3d::edge_t*);

    // The main processor contains the physical qubits and some couplings. The other
    // processor layers contains couplings enabled by TSVs. Each of the processor layers 
    // is guaranteed to be planar. Note that we only need to track planarity for the main
    // processor because the secondary layers only have vertices that are TSVs, and each
    // of these can only have one connection.
    graph::CouplingGraph    get_main_processor(void) { return main_processor; }
    graph::CouplingGraph    get_processor_layer(uint k) { return processor_layers[k]; }
    uint                    get_thickness(void) { return processor_layers.size(); }

    void    compute_coupling_lengths(void);

    std::map<proc3d::edge_t*, fp_t> get_coupling_lengths(void) {
        update_state();
        return coupling_lengths; 
    }

    fp_t get_mean_coupling_length() {
        fp_t sum = 0.0;
        for (auto pair : coupling_lengths) {
            sum += pair.second;   
        }
        return sum / coupling_lengths.size();
    }
private:
    int     find_open_layer(proc3d::vertex_t*, proc3d::vertex_t*, std::set<uint> exclude);  
                                                    // Gets open layer for edge. Returns -1
                                                    // if there are no processor layers 
                                                    // available.
    bool    add_edge_to_processor_layer(proc3d::edge_t*, uint);
    void    remove_coupling_definition(proc3d::edge_t*);
                                                    // Removes coupling definition in
                                                    // coupling_to_edges.

    graph::CouplingGraph                main_processor;
    std::vector<graph::CouplingGraph>   processor_layers;
    uint                                max_thickness;

    std::map<proc3d::vertex_t*, std::vector<proc3d::vertex_t*>> vertex_to_tsv_junctions;   
                                                    // Maps a vertex on the main processor
                                                    // to the TSV vertices above it on the
                                                    // processor layers.
    std::map<proc3d::edge_t*, std::vector<proc3d::edge_t*>>  coupling_to_edges;
                                                    // Maps each coupling of the processor
                                                    // to the set of edges required to
                                                    // implement it.
                                                    // 
                                                    // For example, a coupling that is 
                                                    // on the main processor is mapped to
                                                    // itself (only one edge).
                                                    //
                                                    // A coupling requiring TSVs will
                                                    // map three edges (two vertical,
                                                    // one horizontal).
    std::map<proc3d::edge_t*, fp_t> coupling_lengths;
};

}   // protean
}   // qontra

#endif  // PROTEAN_PROC3D_h
