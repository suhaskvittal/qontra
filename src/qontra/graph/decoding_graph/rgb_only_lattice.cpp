/*
 *  author: Suhas Vittal
 *  date:   27 May 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/decoder/concat_mwpm.h"

namespace qontra {
namespace graph {

using namespace decoding;

sptr<vertex_t>
make_edge_vertex(
        DecodingGraph& gr, sptr<vertex_t> x, sptr<vertex_t> y, evmap_t& evm, uint64_t& idctr) 
{
    auto xy = make_ev_pair(x,y);
    if (!evm.count(xy)) {
        sptr<vertex_t> v = gr.make_and_add_vertex(idctr++);
        evm.put(xy, v);
        return v;
    } else {
        return evm.at(xy);
    }
}

sptr<hyperedge_t>
translate_edge_to_lattice(
        DecodingGraph& gr, sptr<hyperedge_t> e, int color, evmap_t& evm, uint64_t& idctr) 
{
    if (e->get_order() == 2) {
        // Measurement and flag edges.
        auto x = e->get<vertex_t>(0),
             y = e->get<vertex_t>(1);
        if (x->color != color || y->color != color) {
            return nullptr;
        }
        if (!gr.contains(x)) gr.add_vertex(x);
        if (!gr.contains(y)) gr.add_vertex(y);
        // We can just use the edge itself.
        return e;
    } else if (e->get_order() == 3) {
        // Data qubit errors.
        sptr<vertex_t> v;
        std::vector<sptr<vertex_t>> offc;
        for (sptr<vertex_t> x : e->get<vertex_t>()) {
            if (x->color == color) v = x;
            else offc.push_back(x);
        }
        auto x = offc[0],
             y = offc[1];
        if (v->is_boundary_vertex && x->is_boundary_vertex && y->is_boundary_vertex) {
            std::cout << "Found boundary edge.\n";
        }
        if (v->id == 9) {
            std::cout << "Connecting 9(c=0) to " << print_v(x) << ", " << print_v(y) << "\n";
        }
        if (!gr.contains(v)) gr.add_vertex(v);
        auto w = make_edge_vertex(gr, x, y, evm, idctr);
        sptr<hyperedge_t> _e = gr.make_edge({v, w});
        _e->probability = e->probability;
        _e->power = e->power;
        _e->flags = e->flags;
        _e->frames = e->frames;
        return _e;
    } else if (e->get_order() == 4) {
        // Flag edge between two edges of a given color. If it does not follow this format,
        // ignore it. Pre-condition: all weight-4 edges should be saved for last.
        sptr<vertex_t> a1 = nullptr, a2 = nullptr, b1 = nullptr, b2 = nullptr;
        std::vector<int> compl_colors = get_complementary_colors_to({color}, 3);
        int c1 = compl_colors[0],
            c2 = compl_colors[1];
        for (sptr<vertex_t> x : e->get<vertex_t>()) {
            if (x->color == c1) {
                if (a1 == nullptr) a1 = x;
                else                a2 = x;
            } else if (x->color == c2) {
                if (b1 == nullptr) b1 = x;
                else                b2 = x;
            }
        }
        if (a1 == nullptr || a2 == nullptr || b1 == nullptr || b2 == nullptr) {
            return nullptr;
        }
        sptr<vertex_t> v, w;
        if (evm.count(make_ev_pair(a1, b1))) {
            v = evm.at(make_ev_pair(a1,b1));
            w = evm.at(make_ev_pair(a2,b2));
        } else {
            v = evm.at(make_ev_pair(a1,b2));
            w = evm.at(make_ev_pair(a2,b1));
        }
        sptr<hyperedge_t> _e = gr.make_edge({v,w});
        _e->probability = e->probability;
        _e->power = e->power;
        _e->flags = e->flags;
        _e->frames = e->frames;
        return _e;
    } else {
        return nullptr;
    }
}

DecodingGraph
DecodingGraph::make_rgb_only_lattice(int color, evmap_t& evm) {
    DecodingGraph gr;
    // Add edges to rgb only lattice.
    std::vector<sptr<hyperedge_t>> tentative_edges;
    std::vector<sptr<hyperedge_t>> ord4_edges;

    uint64_t idctr = 1L << 48;
    for (sptr<hyperedge_t> e : all_edges) {
        if (e->get_order() == 0) {
            continue;
        } else if (e->get_order() == 4) {
            ord4_edges.push_back(e);
        } else {
            auto _e = translate_edge_to_lattice(gr, e, color, evm, idctr);
            if (_e != nullptr) {
                tentative_edges.push_back(_e);
            }
        }
    }
    // Also handle implicit edge between the boundaries.
    auto ime = get_edge(std::vector<sptr<vertex_t>>{
                    get_boundary_vertex(0),
                    get_boundary_vertex(1),
                    get_boundary_vertex(2)});
    tentative_edges.push_back( translate_edge_to_lattice(gr, ime, color, evm, idctr) );
    for (sptr<hyperedge_t> e : ord4_edges) {
        auto _e = translate_edge_to_lattice(gr, e, color, evm, idctr);
        if (_e != nullptr) {
            tentative_edges.push_back(_e);
        }
    }
    gr.resolve_edges(tentative_edges, 2);
    gr.nod_edges = nod_edges;
    return gr;
}

}   // graph
}   // qontra
