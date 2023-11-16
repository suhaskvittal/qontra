/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 *
 *  Why did I make my own graph class + functions?
 *  Boost has bugs with MPI. So this is stable. I'm
 *  sure I'll regret my decision at some point.
 * */

#ifndef GRAPH_h
#define GRAPH_h

#include "defs.h"

#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

namespace qontra {
namespace graph {

namespace base {

// Templates for the Graph class below
// should subclass base::vertex_t and
// base::edge_t
struct vertex_t {
    uint64_t id;
};

struct edge_t {
    void* src;
    void* dst;
    bool is_undirected=true;
};

}   // base

template <class V_t, class E_t>
class Graph {
public:
    Graph(void)
        :vertices(),
        edges(),
        adjacency_matrix(),
        adjacency_lists(),
        r_adjacency_lists(),
        id_to_vertex(),
        graph_has_changed(false),
        dealloc_on_delete(true)
    {}

    Graph(const Graph& other)
        :vertices(other.vertices), 
        edges(other.edges),
        adjacency_matrix(other.adjacency_matrix),
        adjacency_lists(other.adjacency_lists),
        r_adjacency_lists(other.r_adjacency_lists),
        id_to_vertex(other.id_to_vertex),
        graph_has_changed(other.graph_has_changed),
        dealloc_on_delete(other.dealloc_on_delete)
    {}
    
    virtual ~Graph(void) {
        if (dealloc_on_delete) {
            for (auto v : vertices) delete v;
            for (auto e : edges)    delete e;
        }
    }

    virtual void
    change_id(V_t* v, uint64_t to) {
        id_to_vertex.erase(v->id);
        v->id = to;
        id_to_vertex[to] = v;
    }

    virtual bool
    contains(uint64_t id) {
        return id_to_vertex.count(id);
    }

    virtual bool
    contains(V_t* v) {              // O(1) operation
        return id_to_vertex.count(v->id);
    }

    virtual bool
    contains(E_t* e) {
        auto v1 = (V_t*)e->src;
        auto v2 = (V_t*)e->dst;
        return adjacency_matrix.count(v1)
                && adjacency_matrix[v1].count(v2)
                && adjacency_matrix[v1][v2] == e;
    }

    virtual bool
    contains(V_t* v1, V_t* v2) {    // O(1) operation
        return adjacency_matrix.count(v1)
                && adjacency_matrix[v1].count(v2) 
                && (adjacency_matrix[v1][v2] != nullptr);
    }

    virtual bool
    add_vertex(V_t* v) {            // O(1) operation
        if (contains(v) || id_to_vertex.count(v->id))   return false;
        id_to_vertex[v->id] = v;
        vertices.push_back(v);
        graph_has_changed = true;
        return true;
    }

    virtual bool
    add_edge(E_t* e) {              // O(1) operation
        auto src = (V_t*)e->src;
        auto dst = (V_t*)e->dst;
        if (src == dst)                         return false;
        if (!contains(src) || !contains(dst))   return false;
        if (contains(src, dst))                 return false;
        edges.push_back(e);

        tlm::put(adjacency_matrix, src, dst, e);
        if (e->is_undirected) {
            tlm::put(adjacency_matrix, dst, src, e);
        }

        adjacency_lists[src].push_back(dst);
        r_adjacency_lists[dst].push_back(src);
        if (e->is_undirected) {
            adjacency_lists[dst].push_back(src);
            r_adjacency_lists[src].push_back(dst);
        }

        graph_has_changed = true;
        return true;
    }

    virtual V_t*
    get_vertex(uint64_t id) {        // O(1) operation
        if (!id_to_vertex.count(id))    return nullptr;
        return id_to_vertex[id];
    }

    virtual E_t*
    get_edge(V_t* v1, V_t* v2) {    // O(1) operation
        if (!adjacency_matrix[v1].count(v2))    return nullptr;
        return adjacency_matrix[v1][v2];
    }

    virtual void
    delete_vertex(V_t* v) {         // O(n) operation
        if (!contains(v))   return;
        for (auto it = vertices.begin(); it != vertices.end();) {
            if (*it == v)   it = vertices.erase(it);
            else            it++;
        }

        for (auto it = edges.begin(); it != edges.end();) {
            auto e = *it;
            auto u1 = (V_t*)e->src;
            auto u2 = (V_t*)e->dst;
            if (u1 == v || u2 == v) { 
                // Delete v from the adjacency list of the other vertex.
                V_t* other = u1;
                if (u1 == v)    other = u2;

                auto& adj = adjacency_lists[other];
                for (auto it = adj.begin(); it != adj.end();) {
                    if (*it == v)   it = adj.erase(it);
                    else            it++;
                }

                auto& r_adj = r_adjacency_lists[other];
                for (auto it = r_adj.begin(); it != r_adj.end();) {
                    if (*it == v)   it = adj.erase(it);
                    else            it++;
                }

                adjacency_matrix[u1][u2] = nullptr;
                if (e->is_undirected) {
                    adjacency_matrix[u2][u1] = nullptr;
                }
                // Now delete the edge itself.
                it = edges.erase(it); 
                if (dealloc_on_delete)  delete e;
            } else {
                it++;
            }
        }
        adjacency_lists.erase(v);
        r_adjacency_lists.erase(v);
        graph_has_changed = true;
        if (dealloc_on_delete)  delete v;
    }

    virtual void
    delete_edge(E_t* e) {  // O(m) operation
        if (!contains(e))   return;
        auto src = (V_t*)e->src;
        auto dst = (V_t*)e->dst;

        adjacency_matrix[src][dst] = nullptr;
        
        auto& adj_src = adjacency_lists[src];
        for (auto it = adj_src.begin(); it != adj_src.end();) {
            if (*it == dst) it = adj_src.erase(it);
            else            it++;
        }

        auto& r_adj_dst = r_adjacency_lists[dst];
        for (auto it = r_adj_dst.begin(); it != r_adj_dst.end();) {
            if (*it == src) it = r_adj_src.erase(it);
            else            it++;
        }

        if (e->is_undirected) {
            adjacency_matrix[dst][src] = nullptr;
            auto& adj_dst = adjacency_lists[dst];
            for (auto it = adj_dst.begin(); it != adj_dst.end();) {
                if (*it == src) it = adj_dst.erase(it);
                else            it++;
            }

            auto& r_adj_src = r_adjacency_lists[src];
            for (auto it = r_adj_dst.begin(); it != r_adj_dst.end();) {
                if (*it == dst) it = r_adj_src.erase(it);
                else            it++;
            }
        }

        for (auto it = edges.begin(); it != edges.end();) {
            if (*it == e)   it = edges.erase(it);
            else            it++;
        }

        graph_has_changed=true;
        if (dealloc_on_delete)  delete e;
    }

    std::vector<V_t*>
    get_common_neighbors(V_t* v, V_t* w) {
        auto v_adj = get_neighbors(v);
        auto w_adj = get_neighbors(w);
        for (auto it = v_adj.begin(); it != v_adj.end(); ) {
            if (std::find(w_adj.begin(), w_adj.end(), *it) == w_adj.end()) {
                it = v_adj.erase(it);
            } else {
                it++;
            }
        }
        return v_adj;
    }

    std::vector<V_t*>   get_vertices(void) { return vertices; }
    std::vector<E_t*>   get_edges(void) { return edges; }

    std::vector<V_t*>   get_neighbors(V_t* v) { return adjacency_lists[v]; }
    std::vector<V_t*>   get_incoming(V_t* v) { return r_adjacency_lists[v]; }
    std::vector<V_t*>   get_outgoing(V_t* v) { return adjacency_lists[v]; }

    uint                get_degree(V_t* v) { return get_neighbors(v).size(); }
    uint                get_indegree(V_t* v) { return get_incoming(v).size(); }
    uint                get_outdegree(V_t* v) { return get_degree(v); }
    uint                get_inoutdegree(V_t* v) { return get_indegree(v) + get_outdegree(v); }

    fp_t                get_mean_connectivity(void)
                            { return 2 * ((fp_t)edges.size()) / ((fp_t)vertices.size()); }
    uint                get_max_connectivity(void)
                            { update_state(); return max_degree; }

    bool    dealloc_on_delete;  // Deletes vertices and edges on delete functions if set.

    void    force_update_state(void) { graph_has_changed = true; update_state(); }
protected:
    // Updates graph if graph_has_changed is set.
    // Subclasses should override this method if they track state in any way.
    // 
    // Returns true if the state was updated.
    virtual bool update_state(void) {
        if (!graph_has_changed) return false;
        max_degree = 0;
        for (auto v : vertices) {
            if (get_degree(v) > max_degree) max_degree = get_degree(v);
        }
        graph_has_changed = false;
        return true;
    }

    std::vector<V_t*>   vertices;
    std::vector<E_t*>   edges;

    TwoLevelMap<V_t*, V_t*, E_t*>       adjacency_matrix;
    std::map<V_t*, std::vector<V_t*>>   adjacency_lists;

    // For directed graphs, maintain reverse adjacency lists as well.
    std::map<V_t*, std::vector<V_t*>>   r_adjacency_lists;

    std::map<uint64_t, V_t*>            id_to_vertex;

    bool    graph_has_changed;  // Tracks if the graph has changed. May be useful
                                // for subclasses that need to track state.
private:
    uint max_degree;
};

// Evaluation functions:

template <class V_t>
using ewf_t = std::function<fp_t(V_t*, V_t*)>;

template <class V_t> inline ewf_t<V_t>
unit_ewf_t(void) { return ([] (V_t* v1, V_t* v2) { return 1.0; }); }

}   // graph
}   // qontra

#endif  // GRAPH_h
