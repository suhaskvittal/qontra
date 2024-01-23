/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 *
 *  Why did I make my own graph class + functions?
 *  Boost has bugs with MPI. So this is stable. I'm
 *  sure I'll regret my decision at some point.
 * */

#ifndef QONTRA_GRAPH_h
#define QONTRA_GRAPH_h

#include "qontra/defs.h"

#include <vtils/two_level_map.h>

#include <map>

namespace qontra {
namespace graph {

namespace base {

// Templates for the Graph class below
// should subclass base::vertex_t and
// base::edge_t.

struct vertex_t {
    uint64_t id;
};

struct edge_t {
    sptr<void> src;
    sptr<void> dst;
    bool is_undirected=true;

    template <class V=void> sptr<V> get_source(void);
    template <class V=void> sptr<V> get_target(void);
};

}   // base

// Printing functions for vertices and edges. The user should
// specialize the functions for their own classes.

template <class V> std::string                  print_v(sptr<V>);
template <> std::string                         print_v(sptr<void>);
template <class V=void, class E> std::string    print_e(sptr<E>);

template <class V, class E>
class Graph {
public:
    Graph(void)
        :vertices(),
        edges(),
        adjacency_matrix(),
        adjacency_lists(),
        r_adjacency_lists(),
        id_to_vertex(),
        graph_has_changed(false)
    {}

    Graph(const Graph& other)
        :vertices(other.vertices), 
        edges(other.edges),
        adjacency_matrix(other.adjacency_matrix),
        adjacency_lists(other.adjacency_lists),
        r_adjacency_lists(other.r_adjacency_lists),
        id_to_vertex(other.id_to_vertex),
        graph_has_changed(other.graph_has_changed)
    {}

    virtual void    change_id(sptr<V>, uint64_t to);
    virtual void    manual_update_id(sptr<V>, uint64_t old_id, uint64_t new_id);

    virtual bool    contains(uint64_t id);
    virtual bool    contains(sptr<V>);
    virtual bool    contains(sptr<V>, sptr<V>);
    virtual bool    contains(sptr<E>);

    virtual sptr<V> make_vertex(void);
    virtual sptr<E> make_edge(sptr<V>, sptr<V>, bool is_undirected=true);

    virtual bool    add_vertex(sptr<V>);
    virtual bool    add_edge(sptr<E>);

    virtual sptr<V> make_and_add_vertex(uint64_t id);
    virtual sptr<E> make_and_add_edge(sptr<V>, sptr<V>, bool is_undirected=true);

    virtual sptr<V> get_vertex(uint64_t);
    virtual sptr<E> get_edge(sptr<V>, sptr<V>);
    virtual sptr<E> get_edge(uint64_t, uint64_t);

    virtual void    delete_vertex(sptr<V>);
    virtual void    delete_edge(sptr<E>);

    std::vector<sptr<V>>    get_vertices(void);
    std::vector<sptr<E>>    get_edges(void);

    size_t  n(void);
    size_t  m(void);

    std::vector<sptr<V>>    get_neighbors(sptr<V>);
    std::vector<sptr<V>>    get_incoming(sptr<V>);
    std::vector<sptr<V>>    get_outgoing(sptr<V>);

    std::vector<sptr<V>>    get_common_neighbors(sptr<V>, sptr<V>);

    size_t  get_degree(sptr<V>);
    size_t  get_indegree(sptr<V>);
    size_t  get_outdegree(sptr<V>);
    size_t  get_inoutdegree(sptr<V>);

    fp_t    get_mean_degree(void);
    size_t  get_max_degree(void);

    void    force_update_state(void);
protected:
    // Updates graph if graph_has_changed is set.
    // Subclasses should override this method if they track state in any way.
    // 
    // Returns true if the state was updated.
    virtual bool    update_state(void);

    std::vector<sptr<V>>   vertices;
    std::vector<sptr<E>>   edges;

    vtils::TwoLevelMap<sptr<V>, sptr<V>, sptr<E>>   adjacency_matrix;
    std::map<sptr<V>, std::vector<sptr<V>>>         adjacency_lists;
    // For directed graphs, maintain reverse adjacency lists as well.
    std::map<sptr<V>, std::vector<sptr<V>>>         r_adjacency_lists;

    std::map<uint64_t, sptr<V>>   id_to_vertex;

    bool    graph_has_changed;  // Tracks if the graph has changed. May be useful
                                // for subclasses that need to track state.
private:
    size_t max_degree;
};

}   // graph
}   // qontra

#include "graph.inl"

#endif  // QONTRA_GRAPH_h
