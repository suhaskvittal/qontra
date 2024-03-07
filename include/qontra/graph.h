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

#include <vtils/bijective_map.h>
#include <vtils/two_level_map.h>

#include <map>
#include <vector>

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
    Graph(void) = default;
    Graph(Graph&& other) = default;

    virtual void    change_id(sptr<V>, uint64_t to);
    virtual void    manual_update_id(sptr<V>, uint64_t old_id, uint64_t new_id);

    virtual bool    contains(uint64_t id) const;
    virtual bool    contains(sptr<V>) const;
    virtual bool    contains(sptr<V>, sptr<V>) const;
    virtual bool    contains(sptr<E>) const;

    virtual sptr<V> make_vertex(void) const;
    virtual sptr<V> make_vertex(uint64_t) const;
    virtual sptr<E> make_edge(sptr<V>, sptr<V>, bool is_undirected=true) const;

    virtual bool    add_vertex(sptr<V>);
    virtual bool    add_edge(sptr<E>);

    virtual sptr<V> make_and_add_vertex(uint64_t id);
    virtual sptr<E> make_and_add_edge(sptr<V>, sptr<V>, bool is_undirected=true);

    virtual sptr<V> get_vertex(uint64_t) const;
    virtual sptr<E> get_edge(sptr<V>, sptr<V>) const;
    virtual sptr<E> get_edge(uint64_t, uint64_t) const;

    virtual void    delete_vertex(sptr<V>);
    virtual void    delete_edge(sptr<E>);

    std::vector<sptr<V>>    get_vertices(void) const;
    std::vector<sptr<E>>    get_edges(void) const;

    const vtils::BijectiveMap<sptr<V>, size_t>& get_enumeration_map(void) const;

    size_t  n(void) const;
    size_t  m(void) const;

    std::vector<sptr<V>>    get_neighbors(sptr<V>) const;
    std::vector<sptr<V>>    get_incoming(sptr<V>) const;
    std::vector<sptr<V>>    get_outgoing(sptr<V>) const;

    std::vector<sptr<V>> get_common_neighbors(std::vector<sptr<V>>) const;

    size_t  get_degree(sptr<V>) const;
    size_t  get_indegree(sptr<V>) const;
    size_t  get_outdegree(sptr<V>) const;
    size_t  get_inoutdegree(sptr<V>) const;

    fp_t    get_mean_degree(void) const;
    size_t  get_max_degree(void);

    size_t  const_get_max_degree(void) const;

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

    std::map<uint64_t, sptr<V>> id_to_vertex;
    vtils::BijectiveMap<sptr<V>, size_t> vertex_enum_map;

    bool    graph_has_changed;  // Tracks if the graph has changed. May be useful
                                // for subclasses that need to track state.
private:
    size_t max_degree;
};

}   // graph
}   // qontra

#include "graph.inl"

#endif  // QONTRA_GRAPH_h
