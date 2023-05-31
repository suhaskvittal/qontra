/*
 *  author: Suhas Vittal
 *  date:   31 May 2023
 * */

#include "protean/proc3d.h"

namespace qontra {
namespace protean {

bool
Processor3D::add_edge(proc3d::edge_t* e, bool is_undirected) {
    if (!CouplingGraph::add_edge(e, is_undirected)) return false;
    // Add the edge to the processor while preserving planarity.
    bool will_be_planar = main_processor.test_planarity_after_add(e);
    if (will_be_planar) {
        main_processor.add_edge(e, true);
    } else {
        std::set<uint> open_layers;
        for (uint i = 0; i < get_thickness(); i++)  open_layers.insert(i);
        for (auto tsv : vertex_to_tsv_junctions(e->src)) {
            open_layers.erase(tsv->processor_layer);
        }
        for (auto tsv : vertex_to_tsv_junctions(e->dst)) {
            open_layers.erase(tsv->processor_layer);
        }
        // If there are no open layers, add a new processor layer.
        CouplingMap chip;
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
        proc3d::vertex_t junc1 = new proc3d::vertex_t;
        proc3d::vertex_t junc2 = new proc3d::vertex_t;
        junc1->id = e->src->id | (layer << 28); // Reserve last four bits to identify the layer.
        junc1->processor_layer = layer;
        junc2 = e->dst->id | (layer << 28);
        junc2->processor_layer = layer;
         
        proc3d::edge_t junc1_junc2 = new proc3d::edge_t;
        junc1_junc2->src = junc1;
        junc1_junc2->dst = junc2;
        junc1_junc2->is_vertical = false;

        proc3d::edge_t src_junc1 = new proc3d::edge_t;
        src_junc1->src = e->src;
        src_junc1->dst = junc1;
        src_junc1->is_vertical = true;

        proc3d::edge_t junc2_dst = new proc3d::edge_t;
        junc2_dst->src = junc2;
        junc2_dst->dst = e->ddst;
        junc2_dst->is_vertical = true;

        add_vertex(junc1);
        add_vertex(junc2);
        Graph::add_edge(src_junc1);
        Graph::add_edge(junc1_junc2);
        Graph::add_edge(junc2_dst);

        chip.add_vertex(junc1);
        chip.add_vertex(junc2);
        chip.add_edge(junc1_junc2);

        vertex_to_tsv_junctions[e->src].push_back(junc1);
        vertex_to_tsv_junctions[e->dst].push_back(junc2);
    }
    return true;
}

}   // protean
}   // qontra
