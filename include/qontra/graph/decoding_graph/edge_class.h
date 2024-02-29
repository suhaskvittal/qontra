/*
 *  author: Suhas Vittal
 *  date:   28 February 2024
 * */

#ifndef QONTRA_DECODING_GRAPH_EDGE_CLASS_h
#define QONTRA_DECODING_GRAPH_EDGE_CLASS_h

#include "qontra/graph/decoding_graph.h"

namespace qontra {
namespace graph {

class EdgeClass {
public:
    std::vector<EdgeClass>
        from_edges(const std::vector<sptr<decoding::hyperedge_t>>&);

    void add_vertex(sptr<decoding::vertex_t>);

    sptr<decoding::hyperedge_t> get_representative(void) const;
    std::vector<sptr<decoding::hyperedge_t>> get_edges(void) const;
private:
    EdgeClass(sptr<decoding::hyperedge_t>, std::vector<sptr<decoding::hyperedge_t>>);

    sptr<decoding::hyperedge_t> rep;
    std::vector<sptr<decoding::hyperedge_t>> edges;
};

bool are_in_same_class(sptr<decoding::hyperedge_t>, sptr<decoding::hyperedge_t>);

}   // graph
}   // qontra

#include "edge_class.inl"

#endif  // QONTRA_DECODING_GRAPH_EDGE_CLASS_h
