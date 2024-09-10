/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/edge_class.h"

#include <initializer_list>

#include <vtils/set_algebra.h>

namespace qontra {
namespace graph {

using namespace decoding;

std::vector<sptr<vertex_t>>
DecodingGraph::get_complementary_boundaries_to(std::vector<sptr<vertex_t>> vlist) {
    std::vector<int> colors_in_vlist;
    for (sptr<vertex_t> v : vlist) colors_in_vlist.push_back(v->color);
    std::vector<sptr<vertex_t>> boundary_list;
    for (int c : get_complementary_colors_to(colors_in_vlist, number_of_colors)) {
        boundary_list.push_back(get_boundary_vertex(c));
    }
    return boundary_list;
}

sptr<hyperedge_t>
DecodingGraph::get_best_shared_edge(std::vector<sptr<vertex_t>> vlist) {
    auto common = get_common_hyperedges(vlist);
    if (common.empty()) return nullptr;
    if (common.size() == 1) return common[0];
    // Choose hyperedge with lowest order. Break ties with probability.
    auto it = std::min_element(common.begin(), common.end(),
                [] (sptr<hyperedge_t> x, sptr<hyperedge_t> y)
                {
                    return (x->get_order() < y->get_order())
                        || (x->get_order() == y->get_order() && x->probability > y->probability);
                });
    return *it;
}

sptr<hyperedge_t>
DecodingGraph::get_base_edge(sptr<hyperedge_t> e) {
    if (!edge_class_map.count(e)) {
        return e;
    }
    // We want to (1) get the representative of e's class, (2) flatten the rep,
    // and (3) get the class of the flattened rep. Then, we will search the class
    // for an edge with identical frames to e and an identical number of flags to e.
    const EdgeClass& c = edge_class_map.at(e);
    sptr<hyperedge_t> rep = c.get_representative();
    // Cannot flatten an edge represented by a flag edge.
    if (rep->flags.size()) {
        return e;
    }

    std::vector<sptr<vertex_t>> vlist;
    for (size_t i = 0; i < rep->get_order(); i++) {
        vlist.push_back(rep->get<vertex_t>(i)->get_base());
    }
    sptr<hyperedge_t> frep = get_edge(vlist);
    const EdgeClass& fc = edge_class_map.at(frep);
    for (sptr<hyperedge_t> fe : fc.get_edges()) {
        if (fe->flags.size() == e->flags.size() && fe->frames == e->frames) {
            return fe;
        }
    }
    // Could not flatten this.
    return e;
}

std::vector<sptr<hyperedge_t>>
DecodingGraph::get_flag_edges() {
    // Select one flag edge from each equivalence class that has the most in common
    // with the active flags.
    std::vector<sptr<hyperedge_t>> edge_list;
    for (EdgeClass& c : edge_classes) {
        sptr<hyperedge_t> best_edge = get_best_flag_edge(c.get_edges());
        if (best_edge != nullptr) {
            edge_list.push_back(best_edge);
        }
    }
    return edge_list;
}

std::unordered_map<sptr<hyperedge_t>, sptr<hyperedge_t>>
DecodingGraph::get_best_rep_map() {
    std::unordered_map<sptr<hyperedge_t>, sptr<hyperedge_t>> rep_map;
    for (uint64_t f : active_flags) {
        for (const EdgeClass& c : flag_class_map[f]) {
            sptr<hyperedge_t> best_e = get_best_flag_edge(c.get_edges());
            if (best_e == nullptr) {
                best_e = c.get_representative();
            }
            rep_map[c.get_representative()] = best_e;
        }
    }
    return rep_map;
}

sptr<hyperedge_t>
DecodingGraph::get_best_flag_edge(std::vector<sptr<hyperedge_t>> edge_list) {
    sptr<hyperedge_t> best_edge = nullptr;
    size_t best_nf = 0;
    fp_t best_p = 0.0;

    for (sptr<hyperedge_t> e : edge_list) {
        std::unordered_set<uint64_t> flag_intersect = vtils::immd_set_intersect(e->flags, active_flags);
        const size_t nf = flag_intersect.size();
        const fp_t p = e->probability;
        if (nf == e->flags.size() && 
            (nf > best_nf || (nf > 0 && nf == best_nf && p > best_p)))
        {
            best_edge = e;
            best_nf = nf;
            best_p = p;
        }
    }
    return best_edge;
}

void
DecodingGraph::build_error_polynomial() {
    sptr<hyperedge_t> e0 = edges[0];
    poly_t pX{1 - e0->probability, e0->probability};
    
    fp_t expectation = 0.0;
    for (size_t i = 1; i < edges.size(); i++) {
        sptr<hyperedge_t> e = edges[i];
        poly_t a(pX.size()+1, 0),
               b(pX.size()+1, 0);
        for (size_t j = 0; j < pX.size(); j++) {
            a[j] = pX[j] * (1 - e->probability);
            b[j] = pX[j] * e->probability;
        }
        pX = std::move(a);
        for (size_t j = 0; j < pX.size(); j++) {
            pX[j] += b[j];
            if (i == edges.size()-1) {
                expectation += j * pX[j];
            }
        }
    }
    error_polynomial = std::move(pX);
    expected_errors = expectation;
}

}   // graph
}   // qontra
