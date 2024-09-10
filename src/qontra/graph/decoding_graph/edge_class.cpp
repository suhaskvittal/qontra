/*
 *  author: Suhas Vittal
 *  date:   28 February 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/edge_class.h"

namespace qontra {
namespace graph {

using namespace decoding;

EdgeClass::EdgeClass(sptr<hyperedge_t> v)
    :rep(v),
    edges({v})
{}

std::vector<EdgeClass>
EdgeClass::from_edges(const std::vector<sptr<hyperedge_t>>& edges) {
    std::vector<EdgeClass> eqs;

    for (sptr<hyperedge_t> e : edges) {
        // Go through the list of equivalence classes and find the corresponding one. If
        // none exist for e, make a new one.
        bool found_eq_class = false;
        for (EdgeClass& c : eqs) {
            if (are_in_same_class(e, c.get_representative())) {
                found_eq_class = true;
                c.edges.push_back(e);
            }
        }
        if (found_eq_class) continue;
        // Make a new equivalence class.
        eqs.emplace_back(e);
    }
    // Now, go back through each equivalence class and update the representatives.
    for (EdgeClass& c : eqs) {
        sptr<hyperedge_t> r = *std::min_element(c.edges.begin(), c.edges.end(),
                                [] (sptr<hyperedge_t> x, sptr<hyperedge_t> y)
                                {
                                    return x->flags.size() < y->flags.size();
                                });
        c.rep = r;
    }
    return eqs;
}

// DecodingGraph functions:

sptr<hyperedge_t>
DecodingGraph::get_best_edge_from_class_of(sptr<hyperedge_t> e) {
    if (!edge_class_map.count(e)) {
        return e;
    }
    const EdgeClass& c = edge_class_map.at(e);
    sptr<hyperedge_t> _e = get_best_flag_edge(c.get_edges());
    return _e == nullptr ? e : _e;
}

EdgeClass
DecodingGraph::get_edge_class(sptr<hyperedge_t> e) {
    return edge_class_map.at(e);
}

}   // graph
}   // qontra
