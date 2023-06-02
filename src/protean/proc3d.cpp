/*
 *  author: Suhas Vittal
 *  date:   31 May 2023
 * */

#include "protean/proc3d.h"

namespace qontra {
namespace protean {

using namespace proc3d;

bool
Processor3D::add_edge(proc3d::edge_t* e, bool is_undirected) {
    if (!__Processor3DParent::add_edge(e, is_undirected)) return false;
    auto src = (vertex_t*)e->src;
    auto dst = (vertex_t*)e->dst;
    // Add the edge to the processor while preserving planarity.
    bool will_be_planar = main_processor.test_planarity_after_add(e);
    if (will_be_planar) {
        main_processor.add_edge(e);
    } else {
        std::set<uint> open_layers;
        for (uint i = 0; i < get_thickness(); i++)  open_layers.insert(i);
        for (auto tsv : vertex_to_tsv_junctions[src]) {
            open_layers.erase(tsv->processor_layer);
        }
        for (auto tsv : vertex_to_tsv_junctions[dst]) {
            open_layers.erase(tsv->processor_layer);
        }
        // If there are no open layers, add a new processor layer.
        graph::CouplingGraph chip;
        uint layer;
        if (open_layers.empty()) {
            layer = get_thickness();
            if (layer == max_thickness) return false;   // Cannot add another layer.
            processor_layers.push_back(chip);
        } else {
            layer = *(open_layers.begin());
            chip = processor_layers[layer];
        }
        // Create TSVs and edges.
        proc3d::vertex_t* junc1 = new proc3d::vertex_t;
        proc3d::vertex_t* junc2 = new proc3d::vertex_t;
        junc1->id = src->id | (layer << 28); // Reserve last four bits to identify the layer.
        junc1->processor_layer = layer;
        junc2->id = dst->id | (layer << 28);
        junc2->processor_layer = layer;
         
        proc3d::edge_t* junc1_junc2 = new proc3d::edge_t;
        junc1_junc2->src = junc1;
        junc1_junc2->dst = junc2;
        junc1_junc2->is_vertical = false;
        junc1_junc2->processor_layer = layer;

        proc3d::edge_t* src_junc1 = new proc3d::edge_t;
        src_junc1->src = src;
        src_junc1->dst = junc1;
        src_junc1->is_vertical = true;

        proc3d::edge_t* junc2_dst = new proc3d::edge_t;
        junc2_dst->src = junc2;
        junc2_dst->dst = dst;
        junc2_dst->is_vertical = true;

        Graph::add_edge(e);

        chip.add_vertex(junc1);
        chip.add_vertex(junc2);
        chip.add_edge(junc1_junc2);

        vertex_to_tsv_junctions[src].push_back(junc1);
        vertex_to_tsv_junctions[dst].push_back(junc2);
        coupling_to_edges[e] = std::vector<proc3d::edge_t*>{src_junc1, junc1_junc2, junc2_dst};
    }
    return true;
}

void
Processor3D::delete_edge(proc3d::edge_t* e) {
    __Processor3DParent::delete_edge(e);
    auto impl = coupling_to_edges[e];
    if (impl.size() > 1) {  // Then this is a complex coupling, so remove the corresponding
                            // metadata.
        auto vert1 = impl[0];
        auto horiz = impl[1];
        auto vert2 = impl[2];

        auto& src_tsv_j = vertex_to_tsv_junctions[(vertex_t*)vert1->src];
        for (auto it = src_tsv_j.begin(); it != src_tsv_j.end();) {
            if (*it == (vertex_t*)vert1->dst)   it = src_tsv_j.erase(it);
            else                                it++;
        }
        auto& dst_tsv_j = vertex_to_tsv_junctions[(vertex_t*)vert2->dst];
        for (auto it = dst_tsv_j.begin(); it != dst_tsv_j.end();) {
            if (*it == (vertex_t*)vert2->src)   it = dst_tsv_j.erase(it);
            else                                it++;
        }
        processor_layers[horiz->processor_layer].delete_edge(horiz);
    }
    coupling_to_edges.erase(e);
}

}   // protean
}   // qontra
