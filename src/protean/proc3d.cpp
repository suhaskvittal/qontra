/* 
 *  author: Suhas Vittal
 *  date:   31 May 2023
 * */

#include "protean/proc3d.h"

#define SQR(x) ((x)*(x))

namespace qontra {
namespace protean {

using namespace proc3d;

bool
Processor3D::add_vertex(vertex_t* v) {
    if (!main_processor.add_vertex(v))  return false;
    if (!Graph::add_vertex(v))  return false;
    return true;
}

bool
Processor3D::add_edge(edge_t* e, bool force_out_of_plane) {
    auto src = (vertex_t*)e->src;
    auto dst = (vertex_t*)e->dst;
    if (!force_out_of_plane) {
        if (!__Processor3DParent::add_edge(e)) return false;
        // Add the edge to the processor while preserving planarity.
        bool will_be_planar = main_processor.test_planarity_after_add(e);
        if (will_be_planar) {
            main_processor.add_edge(e, will_be_planar);
            coupling_to_edges[e] = std::vector<edge_t*>{e};
            return true;
        }
    }
    std::set<uint> not_open;
add_try_placement:
    // Find open layer.
    int layer = find_open_layer(src, dst, not_open);
    // If there are no open layers, add a new processor layer.
    if (layer < 0) {
        layer = get_thickness();
        if (layer == max_thickness) return false;   // Cannot add another layer.
        graph::CouplingGraph new_chip;
        new_chip.dealloc_on_delete = false;
        processor_layers.push_back(new_chip);
    }
    if (!add_edge_to_processor_layer(e, layer)) {
        not_open.insert(layer);
        goto add_try_placement;
    }
    return true;
}

void
Processor3D::delete_vertex(vertex_t* v) {
    for (auto w : get_neighbors(v)) {
        auto e = get_edge(v, w);        
        remove_coupling_definition(e);
    } 
    vertex_to_tsv_junctions.erase(v);
    main_processor.delete_vertex(v);
    __Processor3DParent::delete_vertex(v);
}

void
Processor3D::delete_edge(edge_t* e) {
    if (!contains(e))   return;
    if (e == nullptr)   return;
    remove_coupling_definition(e);
    __Processor3DParent::delete_edge(e);
}

void
Processor3D::reallocate_edges() {
    for (auto e : edges) {
        if (main_processor.contains(e))    continue;
        bool p = main_processor.test_planarity_after_add(e);
        // First, try to add to main processor.
        if (p) {
            remove_coupling_definition(e);
            e->processor_layer = 0;
            coupling_to_edges[e] = {e};
            main_processor.add_edge(e, p);
        } else {
            std::set<uint> not_open;
realloc_try_placement:
            auto impl = coupling_to_edges[e];
            uint curr_layer = impl[1]->processor_layer-1;
            // Check if there are any open layers below the current layer.
            int new_layer = find_open_layer(
                                (vertex_t*)e->src, (vertex_t*)e->dst, not_open);
            if (new_layer < 0 || curr_layer < new_layer)  continue;
            remove_coupling_definition(e);
            if (!add_edge_to_processor_layer(e, new_layer)) {
                not_open.insert(new_layer);
                goto realloc_try_placement;
            }
        }
    }
    // Check if any processor layers are empty. If so, remove them.
    for (auto it = processor_layers.begin(); it != processor_layers.end(); ) {
        if (it->get_edges().empty() || it->get_vertices().empty()) {
            it = processor_layers.erase(it); 
        } else {
            it++;
        }
    }
}

void
Processor3D::force_out_of_plane(edge_t* e) {
    if (!contains(e))                   return;
    if (!main_processor.contains(e))    return;
    remove_coupling_definition(e);
    main_processor.delete_edge(e);
    add_edge(e, true);
}

int
Processor3D::find_open_layer(vertex_t* src, vertex_t* dst,
        std::set<uint> exclude) 
{
    std::set<uint> open_layers;
    for (uint i = 0; i < get_thickness(); i++)  open_layers.insert(i);
    for (auto tsv : vertex_to_tsv_junctions[src]) {
        open_layers.erase(tsv->processor_layer-1);
    }
    for (auto tsv : vertex_to_tsv_junctions[dst]) {
        open_layers.erase(tsv->processor_layer-1);
    }
    for (auto x : exclude)  open_layers.erase(x);
    
    uint layer = open_layers.empty() ? -1 : *open_layers.begin();
    return layer;
}

bool
Processor3D::add_edge_to_processor_layer(edge_t* e, uint layer) {
    auto src = (vertex_t*)e->src;
    auto dst = (vertex_t*)e->dst;
    graph::CouplingGraph& chip = processor_layers[layer];

    layer++;    // Increment so junction id is clearly shows the layer
    // Create TSVs and edges.
    vertex_t* junc1 = new vertex_t;
    vertex_t* junc2 = new vertex_t;
    junc1->id = src->id | (layer << 28); // Reserve last four bits to identify the layer.
    junc1->processor_layer = layer;
    junc2->id = dst->id | (layer << 28);
    junc2->processor_layer = layer;
     
    edge_t* junc1_junc2 = new edge_t;
    junc1_junc2->src = junc1;
    junc1_junc2->dst = junc2;
    junc1_junc2->is_vertical = false;
    junc1_junc2->processor_layer = layer;

    edge_t* src_junc1 = new edge_t;
    src_junc1->src = src;
    src_junc1->dst = junc1;
    src_junc1->is_vertical = true;

    edge_t* junc2_dst = new edge_t;
    junc2_dst->src = junc2;
    junc2_dst->dst = dst;
    junc2_dst->is_vertical = true;

    chip.add_vertex(junc1);
    chip.add_vertex(junc2);
    if (!chip.add_edge(junc1_junc2)) {
        chip.delete_vertex(junc1);
        chip.delete_vertex(junc2);
        delete src_junc1;
        delete junc1_junc2;
        delete junc2_dst;
        return false;
    }

    vertex_to_tsv_junctions[src].push_back(junc1);
    vertex_to_tsv_junctions[dst].push_back(junc2);
    coupling_to_edges[e] = std::vector<edge_t*>{src_junc1, junc1_junc2, junc2_dst};
    return true;
}

void
Processor3D::remove_coupling_definition(edge_t* e) {
    auto impl = coupling_to_edges[e];
    if (impl.size() == 3) {  // Then this is a complex coupling, so remove the corresponding
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
        processor_layers[horiz->processor_layer-1].delete_vertex((vertex_t*)vert1->dst);
        processor_layers[horiz->processor_layer-1].delete_vertex((vertex_t*)vert2->src);
        processor_layers[horiz->processor_layer-1].delete_edge(horiz);

        if (dealloc_on_delete) {
            delete vert1;
            delete vert2;
        }
    }
    coupling_to_edges.erase(e);
}

void
Processor3D::compute_coupling_lengths() {
    lemon::ListGraph gr;
    std::map<vertex_t*, lemon::ListGraph::Node> lemon_node_map;
    for (auto v : main_processor.get_vertices()) {
        auto x = gr.addNode();
        lemon_node_map[(vertex_t*)v] = x;
    }
    for (auto e : main_processor.get_edges()) {
        auto x = lemon_node_map[(vertex_t*)e->src];
        auto y = lemon_node_map[(vertex_t*)e->dst];
        gr.addEdge(x, y);
    }
    lemon::PlanarEmbedding emb(gr);
    if (!emb.run(false))    return;
    lemon::PlanarDrawing drawing(gr);
    drawing.run(emb.embeddingMap());

    for (auto e : edges) {
        // So there is a vertical component to out of plane edges,
        // but we assume this is negligible.
        auto v = (vertex_t*)e->src;
        auto w = (vertex_t*)e->dst;
        auto locv = drawing[lemon_node_map[v]];
        auto locw = drawing[lemon_node_map[w]];
        
        coupling_lengths[e] = sqrt(SQR(locv.x - locw.x) + SQR(locv.y - locw.y));
    }
}

}   // protean
}   // qontra
