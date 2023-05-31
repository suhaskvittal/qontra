/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#ifndef COUPLING_GRAPH_h
#define COUPLING_GRAPH_h

#include "defs.h"
#include "graph/graph.h"

#include <lemon/list_graph.h>

#include <set>
#include <vector>

namespace qontra {
namespace graph {

namespace coupling {

struct vertex_t : base::vertex_t {
    int32_t logical_owner;  // Number of logical qubit which contains the qubit.
};

struct edge_t : base::edge_t<vertex_t> {
};

}   // coupling

#define __CouplingGraphParent Graph<coupling::vertex_t, coupling::edge_t>

class CouplingGraph : public __CouplingGraphParent {
public:
    CouplingGraph(bool track_planarity=false)
        :Graph(), track_planarity(track_planarity)
    {}

    CouplingGraph(const CouplingGraph& other)
        :Graph(other),
        track_planarity(track_planarity),
        is_planar(is_planar)
    {}

    // The functions test the planarity of the coupling graph after making
    // a modification (adding/deleting a group of edges, or deleting a group
    // of vertices).
    //
    // These functions have single operand versions.

    bool    test_planarity_after_add(const std::vector<coupling::edge_t*>&);
    bool    test_planarity_after_delete(const std::vector<coupling::vertex_t*>&);
    bool    test_planarity_after_delete(const std::vector<coupling::edge_t*>&);

    bool    test_planarity_after_add(coupling::edge_t* e) 
                { std::vector s{e}; return test_planarity_after_add(s); }
    bool    test_planarity_after_delete(coupling::vertex_t* v)
                { std::vector s{v}; return test_planarity_after_delete(s); }
    bool    test_planarity_after_delete(coupling::edge_t*)
                { std::vector s{e}; return test_planarity_after_delete(s); }

    // Graph modification operations that can modify planarity now check planarity.
    // We note that if track_planarity is false, then planarity is not checked.

    bool    add_edge(coupling::edge_t* e) override
                { bool out = Graph::add_edge(e); check_planarity(); return out; }
    void    delete_vertex(coupling::vertex_t* v) override
                { Graph::delete_vertex(v); check_planarity(); }
    void    delete_edge(coupling::edge_t* e) override
                { Graph::delete_edge(e); check_planarity(); }

    // Use the below functions to perform planarity modifying operations in a batch.
    // 
    // These functions can be used to update the planarity all at once instead of
    // updating after every individual operation.

    void    add_edges(const std::vector<edge_t*> edges) 
                { for (auto e : edges) { Graph::add_edge(e); } check_planarity(); }
    void    delete_vertices(const std::vector<edge_t*> vertices)
                { for (auto v : vertices) { Graph::delete_vertex(v); } check_planarity(); }
    void    delete_edges(const std::vector<edge_t*> edges)
                { for (auto e : edges) { Graph::delete_edge(e); } check_planarity(); }

    // Use the below functions to update planarity given a prior test
    // (i.e. one might test planarity after adding an edge and only add the edge if the
    // graph remains planar).

    bool    add_edge(coupling::edge_t* e, bool p)
                { is_planar = p; return Graph::add_edge(e); }
    void    delete_vertex(coupling::vertex_t* v, bool p)
                { Graph::delete_vertex(v); is_planar = p; }
    void    delete_edge(coupling::edge_t* e, bool p)
                { Graph::delete_edge(e); is_planar = p; }

    bool    planar(void) { return is_planar; }
    bool    tracking_planarity(void) { return track_planarity; }

    void toggle_track_planarity(void) { 
        track_planarity = !track_planarity;
        if (track_planarity)    check_planarity();
    }
private:
    void    check_planarity(void);  // Uses LEMON

    bool    is_planar;
    bool    track_planarity;
};

}   // graph
}   // qontra

#endif  // COUPLING_GRAPH_h
