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

    template <class V=void> std::vector<sptr<V>> get(void) const;
    template <class V=void> sptr<V>              get(size_t) const;
    
    size_t  get_order(void) const;
};

}   // base

template <class V=void, class HE>
std::string print_he(sptr<HE>);

typedef Graph<base::vertex_t, base::edge_t> IncidenceGraph;

template <class V, class HE>
class HyperGraph {
public:
    HyperGraph(void)
        :vertices(),
        edges(),
        incidence_graph(std::make_unique<IncidenceGraph>()),
        incidence_object_map(),
        incidence_id_map(),
        adjacency_lists(),
        adjacency_mult_map(),
        id_to_vertex(),
        graph_has_changed(false),
        mean_degree(0.0),
        max_degree(0),
        max_order(0)
    {}

    HyperGraph(HyperGraph&& other) = default;
    
    virtual void    change_id(sptr<V>, uint64_t to);
    virtual void    manual_update_id(sptr<V>, uint64_t old_id, uint64_t new_id);

    virtual bool    contains(uint64_t id) const;
    virtual bool    contains(sptr<V>) const;
    virtual bool    contains(sptr<HE>) const;

    virtual bool    contains_(std::vector<sptr<void>>) const;
    virtual bool    contains(std::vector<sptr<V>>) const;

    virtual sptr<V>     make_vertex(void) const;
    virtual sptr<V>     make_vertex(uint64_t) const;
    virtual sptr<HE>    make_edge(std::vector<sptr<V>>) const;

    virtual bool        add_vertex(sptr<V>);
    virtual bool        add_edge(sptr<HE>);

    virtual sptr<V>     make_and_add_vertex(uint64_t id);
    virtual sptr<HE>    make_and_add_edge(std::vector<sptr<V>>);

    virtual sptr<V>     get_vertex(uint64_t) const;
    virtual sptr<HE>    get_edge_(std::vector<sptr<void>>) const;
    virtual sptr<HE>    get_edge(std::vector<sptr<V>>) const;
    virtual sptr<HE>    get_edge(std::vector<uint64_t>) const;
    
    virtual void    delete_vertex(sptr<V>);
    virtual void    delete_edge(sptr<HE>);

    std::vector<sptr<V>>    get_vertices(void) const;
    std::vector<sptr<HE>>   get_edges(void) const;

    size_t  n(void) const;
    size_t  m(void) const;

    bool    has_endpoint(sptr<HE>, sptr<V>) const;

    bool                    share_hyperedge(std::vector<sptr<V>>) const;
    std::vector<sptr<HE>>   get_common_hyperedges(std::vector<sptr<V>>) const;

    std::vector<sptr<V>>    get_neighbors(sptr<V>) const;
    std::vector<sptr<V>>    get_common_neighbors(std::vector<sptr<V>>) const;

    size_t  get_degree(sptr<V>) const;

    fp_t    get_mean_degree(void);
    size_t  get_max_degree(void);
    size_t  get_max_order(void);

    fp_t    const_get_mean_degree(void) const;
    size_t  const_get_max_degree(void) const;
    size_t  const_get_max_order(void) const;

    void    force_update_state(void);
protected:
    virtual bool    update_state(void);

    template <class PTR> sptr<base::vertex_t> add_incidence_vertex(PTR);
    template <class PTR> sptr<base::vertex_t> get_incidence_vertex(PTR) const;

    sptr<V>     get_vertex_from_incidence_vertex(sptr<base::vertex_t>) const;
    sptr<HE>    get_edge_from_incidence_vertex(sptr<base::vertex_t>) const;

    std::vector<sptr<V>>    vertices;
    std::vector<sptr<HE>>   edges;

    uptr<IncidenceGraph> incidence_graph;
    std::map<sptr<base::vertex_t>, sptr<void>> incidence_object_map;
    std::map<sptr<void>, uint64_t> incidence_id_map;
    std::map<sptr<V>, std::vector<sptr<V>>> adjacency_lists;
    vtils::TwoLevelMap<sptr<V>, sptr<V>, size_t> adjacency_mult_map;

    std::map<uint64_t, sptr<V>> id_to_vertex;

    bool graph_has_changed;
private:
    bool    is_edge_valid(sptr<HE>);
    void    update_adjacency_lists_after_delete(sptr<HE>);

    fp_t   mean_degree;
    size_t max_degree;
    size_t max_order;
};

}   // graph
}   // qontra

#include "inl/hypergraph.inl"

#endif  // QONTRA_HYPERGRAPH_h
