/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#ifndef QONTRA_HYPERGRAPH_h
#define QONTRA_HYPERGRAPH_h

#include "qontra/graph.h"

#include <set>

namespace qontra {
namespace graph {

namespace base {

struct hyperedge_t {
    std::vector<sptr<void>> endpoints;

    template <class V=void> std::vector<sptr<V>> operator()(void) const;
    template <class V=void> sptr<V>              operator[](size_t) const;
    
    size_t  get_order(void) const;
};

}   // base

template <class V=void, class HE> std::string   print_he(sptr<HE>);

template <class V, class HE>
class Hypergraph {
public:
    Hypergraph(void)
        :vertices(),
        edges(),
        adjacency_lists(),
        adjacency_mult_map(),
        id_to_vertex(),
        graph_has_changed(false),
        mean_degree(0.0),
        max_degree(0),
        max_order(0)
    {}

    Hypergraph(const Hypergraph& other)
        :vertices(other.vertices),
        edges(other.edges),
        adjacency_lists(other.adjacency_lists),
        adjacency_mult_map(other.adjacency_mult_map),
        id_to_vertex(other.id_to_vertex),
        graph_has_changed(other.graph_has_changed),
        mean_degree(other.mean_degree),
        max_degree(other.max_degree),
        max_order(other.max_order)
    {}
    
    virtual void    change_id(sptr<V>, uint64_t to);
    virtual void    manual_update_id(sptr<V>, uint64_t old_id, uint64_t new_id);

    virtual bool    contains(uint64_t id) const;
    virtual bool    contains(sptr<V>) const;
    virtual bool    contains(sptr<HE>) const;

    template <class CONTAINER> virtual bool contains(CONTAINER) const;
    template <class ITER> virtual bool      contains(ITER begin, ITER end) const;

    virtual sptr<V>                             make_vertex(void) const;
    virtual sptr<V>                             make_vertex(uint64_t) const;
    template <class CONTAINER> virtual sptr<HE> make_edge(CONTAINER) const;
    template <class ITER> virtual sptr<HE>      make_edge(ITER begin, ITER end) const;

    virtual bool    add_vertex(sptr<V>);
    virtual bool    add_edge(sptr<HE>);

    virtual sptr<V>                             make_and_add_vertex(uint64_t id);
    template <class CONTAINER> virtual sptr<HE> make_and_add_edge(CONTAINER);
    template <class ITER> virtual sptr<HE>      make_and_add_edge(ITER begin, ITER end);

    virtual sptr<V>                             get_vertex(uint64_t) const;
    template <class CONTAINER> virtual sptr<HE> get_edge(CONTAINER) const;
    template <class ITER> virtual sptr<HE>      get_edge(ITER, ITER) const;
    template <class CONTAINER> virtual sptr<HE> get_edge_from_ids(CONTAINER) const;
    template <class ITER> virtual sptr<HE>      get_edge_from_ids(ITER, ITER) const;
    
    virtual void    delete_vertex(sptr<V>);
    virtual void    delete_edge(sptr<HE>);

    std::vector<sptr<V>>    get_vertices(void) const;
    std::vector<sptr<HE>>   get_edges(void) const;

    size_t  n(void) const;
    size_t  m(void) const;

    bool    has_endpoint(sptr<HE>, sptr<V>) const;

    template <class CONTAINER> bool                     share_hyperedge(CONTAINER) const;
    template <class CONTAINER> std::vector<sptr<HE>>    get_common_hyperedges(CONTAINER) const;

    std::vector<sptr<V>>                            get_neighbors(sptr<V>) const;
    template <class CONTAINER> std::vector<sptr<V>> get_common_neighbors(CONTAINER) const;

    size_t  get_degree(sptr<V>) const;

    fp_t    get_mean_degree(void);
    size_t  get_max_degree(void);
    size_t  get_max_order(void);

    fp_t    const_get_mean_degree(void) const;
    size_t  const_get_max_degree(void) const;
    size_t  const_get_max_order(void) const;

    void    force_update_state(void);
protected:
    struct inc_vertex_t : base::vertex_t { sptr<void> obj_p; bool is_vertex; };
    struct inc_edge_t : base::edge_t {};
    typedef Graph<inc_vertex_t, inc_edge_t> IncidenceGraph;

    virtual bool    update_state(void);

    template <class PTR>
    sptr<inc_vertex_t> get_incidence_vertex(PTR) const;

    std::vector<sptr<V>>    vertices;
    std::vector<sptr<HE>>   edges;

    IncidenceGraph incidence_graph;
    std::map<sptr<V>, std::vector<sptr<V>>> adjacency_lists;
    TwoLevelMap<sptr<V>, sptr<V>, size_t>   adjacency_mult_map;

    std::map<uint64_t, sptr<V>> id_to_vertex;

    bool graph_has_changed;
private:
    bool    is_edge_valid(sptr<HE>);

    fp_t   mean_degree;
    size_t max_degree;
    size_t max_order;
};

}   // graph
}   // qontra

#include "hypergraph.inl"

#endif  // QONTRA_HYPERGRAPH_h
