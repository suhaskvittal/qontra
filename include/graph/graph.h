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
#include "defs/two_level_map.h"

#include <iostream>
#include <map>

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
    sptr<void> src;
    sptr<void> dst;
    bool is_undirected=true;
};

}   // base

template <class V> inline std::string 
print_v(sptr<V> v) {
    return std::to_string(v->id);
}

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

    virtual inline void
    change_id(sptr<V> v, uint64_t to) {
        id_to_vertex.erase(v->id);
        v->id = to;
        id_to_vertex[to] = v;
    }

    virtual inline bool
    contains(uint64_t id) {
        return id_to_vertex.count(id);
    }

    virtual inline bool
    contains(sptr<V> v) {              // O(1) operation
        return id_to_vertex.count(v->id);
    }

    virtual inline bool
    contains(sptr<E> e) {
        auto v1 = std::reinterpret_pointer_cast<V>(e->src);
        auto v2 = std::reinterpret_pointer_cast<V>(e->dst);
        return adjacency_matrix.count(v1)
                && adjacency_matrix[v1].count(v2)
                && adjacency_matrix[v1][v2] == e;
    }

    virtual inline bool
    contains(sptr<V> v1, sptr<V> v2) {    // O(1) operation
        return adjacency_matrix.count(v1)
                && adjacency_matrix[v1].count(v2) 
                && (adjacency_matrix[v1][v2] != nullptr);
    }

    virtual inline bool
    add_vertex(sptr<V> v) {            // O(1) operation
        if (contains(v) || id_to_vertex.count(v->id)) return false;
        id_to_vertex[v->id] = v;
        vertices.push_back(v);
        graph_has_changed = true;
        return true;
    }

    virtual inline bool
    add_edge(sptr<E> e) {              // O(1) operation
        auto src = std::reinterpret_pointer_cast<V>(e->src);
        auto dst = std::reinterpret_pointer_cast<V>(e->dst);
        if (src == dst)                         return false;
        if (!contains(src) || !contains(dst))   return false;
        if (contains(src, dst))                 return false;
        edges.push_back(e);

        tlm_put(adjacency_matrix, src, dst, e);
        if (e->is_undirected) {
            tlm_put(adjacency_matrix, dst, src, e);
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

    virtual inline sptr<V>
    get_vertex(uint64_t id) {        // O(1) operation
        if (!id_to_vertex.count(id)) return nullptr;
        return id_to_vertex[id];
    }

    virtual inline sptr<E>
    get_edge(sptr<V> v1, sptr<V> v2) {    // O(1) operation
        if (!adjacency_matrix[v1].count(v2)) return nullptr;
        return adjacency_matrix[v1][v2];
    }

    virtual void
    delete_vertex(sptr<V> v) {         // O(n) operation
        if (!contains(v))   return;
        for (auto it = vertices.begin(); it != vertices.end();) {
            if (*it == v)   it = vertices.erase(it);
            else            it++;
        }

        for (auto it = edges.begin(); it != edges.end();) {
            auto e = *it;
            auto u1 = std::reinterpret_pointer_cast<V>(e->src);
            auto u2 = std::reinterpret_pointer_cast<V>(e->dst);
            if (u1 == v || u2 == v) { 
                // Delete v from the adjacency list of the other vertex.
                sptr<V> other = u1;
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
            } else {
                it++;
            }
        }
        adjacency_lists.erase(v);
        r_adjacency_lists.erase(v);
        graph_has_changed = true;
    }

    virtual void
    delete_edge(sptr<E> e) {  // O(m) operation
        if (!contains(e))   return;
        auto src = std::reinterpret_pointer_cast<V>(e->src);
        auto dst = std::reinterpret_pointer_cast<V>(e->dst);

        adjacency_matrix[src][dst] = nullptr;
        
        auto& adj_src = adjacency_lists[src];
        for (auto it = adj_src.begin(); it != adj_src.end();) {
            if (*it == dst) it = adj_src.erase(it);
            else            it++;
        }

        auto& r_adj_dst = r_adjacency_lists[dst];
        for (auto it = r_adj_dst.begin(); it != r_adj_dst.end();) {
            if (*it == src) it = r_adj_dst.erase(it);
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
    }

    std::vector<sptr<V>>
    get_common_neighbors(sptr<V> v, sptr<V> w) {
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

    inline std::vector<sptr<V>> get_vertices(void) { return vertices; }
    inline std::vector<sptr<E>> get_edges(void) { return edges; }

    inline size_t n(void) { return vertices.size(); }
    inline size_t m(void) { return edges.size(); }

    inline std::vector<sptr<V>> get_neighbors(sptr<V> v) { return adjacency_lists[v]; }
    inline std::vector<sptr<V>> get_incoming(sptr<V> v) { return r_adjacency_lists[v]; }
    inline std::vector<sptr<V>> get_outgoing(sptr<V> v) { return adjacency_lists[v]; }

    inline size_t get_degree(sptr<V> v) { return get_neighbors(v).size(); }
    inline size_t get_indegree(sptr<V> v) { return get_incoming(v).size(); }
    inline size_t get_outdegree(sptr<V> v) { return get_degree(v); }
    inline size_t get_inoutdegree(sptr<V> v) { return get_indegree(v) + get_outdegree(v); }

    inline fp_t get_mean_connectivity(void) {
        return 2 * ((fp_t)edges.size()) / ((fp_t)vertices.size()); 
    }

    inline size_t get_max_connectivity(void) {
        update_state();
        return max_degree;
    }

    inline void force_update_state(void) {
        graph_has_changed = true; 
        update_state(); 
    }

    virtual inline sptr<V> make_vertex(void) {
        return std::make_shared<V>();
    }

    virtual inline sptr<E> make_edge(sptr<V> src, sptr<V> dst) {
        sptr<E> e = std::make_shared<E>();
        e->src = src;
        e->dst = dst;
        return e;
    }

    virtual inline sptr<V> get_source(sptr<E> e) {
        return std::reinterpret_pointer_cast<V>(e->src);
    }

    virtual inline sptr<V> get_target(sptr<E> e) {
        return std::reinterpret_pointer_cast<V>(e->dst);
    }
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

    std::vector<sptr<V>>   vertices;
    std::vector<sptr<E>>   edges;

    TwoLevelMap<sptr<V>, sptr<V>, sptr<E>>  adjacency_matrix;
    std::map<sptr<V>, std::vector<sptr<V>>> adjacency_lists;

    // For directed graphs, maintain reverse adjacency lists as well.
    std::map<sptr<V>, std::vector<sptr<V>>> r_adjacency_lists;

    std::map<uint64_t, sptr<V>>   id_to_vertex;

    bool    graph_has_changed;  // Tracks if the graph has changed. May be useful
                                // for subclasses that need to track state.
private:
    size_t max_degree;
};

}   // graph
}   // qontra

#endif  // GRAPH_h
