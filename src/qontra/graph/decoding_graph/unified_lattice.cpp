/*
 *  author: Suhas Vittal
 *  date:   31 March 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/unified_lattice.h"

namespace qontra {
namespace graph {

using namespace decoding;

void
add_edge_to_ufl(
        DecodingGraph& ufl,
        sptr<hyperedge_t> e,
        int c1,
        int c2, 
        std::vector<sptr<hyperedge_t>>& edge_list,
        std::vector<sptr<vertex_t>>& boundary_paired_list) 
{
    const fp_t p = e->probability;
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<vertex_t> v = e->get<vertex_t>(i);
        if (v->color != c1 && v->color != c2) continue;
        for (size_t j = i+1; j < e->get_order(); j++) {
            sptr<vertex_t> w = e->get<vertex_t>(j);
            if (w->color != c1 && w->color != c2) continue;
            if (v->is_boundary_vertex || w->is_boundary_vertex) {
                if (v->is_boundary_vertex && !w->is_boundary_vertex) {
                    sptr<vertex_t> x = ufl.get_vertex(unified_lattice_id(w->id, c1, c2));
                    boundary_paired_list.push_back(x);
                } else if (!v->is_boundary_vertex && w->is_boundary_vertex) {
                    sptr<vertex_t> x = ufl.get_vertex(unified_lattice_id(v->id, c1, c2));
                    boundary_paired_list.push_back(x);
                }
                continue;
            }
            const uint64_t xid = unified_lattice_id(v->id, c1, c2),
                            yid = unified_lattice_id(w->id, c1, c2);
            sptr<vertex_t> x = ufl.get_vertex(xid),
                            y = ufl.get_vertex(yid);
            sptr<hyperedge_t> f = ufl.make_edge({x, y});
            f->probability = p;
            f->flags = e->flags;
            edge_list.push_back(f);
        }
    }
}

DecodingGraph::DecodingGraph() 
    :number_of_colors(0),
    renorm_factor(1.0)
{}

DecodingGraph
DecodingGraph::make_unified_lattice(
        std::map<sptr<vertex_t>, sptr<vertex_t>>& ufl_map,
        std::map<sptr<hyperedge_t>, sptr<vertex_t>>& crease_edge_map)
{
    DecodingGraph ufl;
    // Copy all non boundary vertices over to ufl.
    for (sptr<vertex_t> v : get_vertices()) {
        if (v->is_boundary_vertex) continue;
        // Each vertex will have two copies.
        for (int c1 = 0; c1 < number_of_colors; c1++) {
            for (int c2 = c1+1; c2 < number_of_colors; c2++) {
                if (v->color != c1 && v->color != c2) continue;
                sptr<vertex_t> x = 
                    ufl.make_and_add_vertex(unified_lattice_id(v->id, c1, c2));
                x->base = x;
                ufl_map[x] = v;
            }
        }
    }
    // Accumulate all tentative edges.
    std::vector<sptr<hyperedge_t>> tentative_edges;
    // First go through normal edges.
    for (sptr<hyperedge_t> e : all_edges) {
        if (e->get_order() == 0 || e->get_order() > 3) continue;
        // Track would-be boundary edges. These will be merged together in the unified lattice.
        std::vector<sptr<vertex_t>> boundary_paired_list;
        for (int c1 = 0; c1 < number_of_colors; c1++) {
            for (int c2 = c1+1; c2 < number_of_colors; c2++) {
                add_edge_to_ufl(ufl, e, c1, c2, tentative_edges, boundary_paired_list);
            }
        }
        // Make crease edges.
        for (size_t i = 0; i < boundary_paired_list.size(); i++) {
            sptr<vertex_t> x = boundary_paired_list.at(i);
            for (size_t j = i+1; j < boundary_paired_list.size(); j++) {
                sptr<vertex_t> y = boundary_paired_list.at(j);
                sptr<hyperedge_t> f = ufl.make_edge({x, y});
                f->probability = e->probability;
                crease_edge_map[f] = boundary_paired_list.size() / 2;
                tentative_edges.push_back(f);
            }
        }
    }
    ufl.resolve_edges(tentative_edges, 2);
    ufl.nod_edges = nod_edges;
    return ufl;
}

}   // graph
}   // qontra
