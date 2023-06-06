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
    uint id;
};

struct edge_t {
    void* src;
    void* dst;
};

}   // base

template <class V_t, class E_t>
class Graph {
public:
    Graph(void)
        :vertices(), edges(), adjacency_matrix(), adjacency_lists(), id_to_vertex()
    {}

    Graph(const Graph& other)
        :vertices(other.vertices), 
        edges(other.edges),
        adjacency_matrix(other.adjacency_matrix),
        adjacency_lists(other.adjacency_lists),
        id_to_vertex(other.id_to_vertex)
    {}
    
    virtual ~Graph(void) {}

    virtual bool
    contains(V_t* v) {              // O(1) operation
        return id_to_vertex.count(v->id);
    }

    virtual bool
    contains(E_t* e) {
        auto v1 = (V_t*)e->src;
        auto v2 = (V_t*)e->dst;
        return adjacency_matrix[v1].count(v2)
                && adjacency_matrix[v1][v1] == e;
    }

    virtual bool
    contains(V_t* v1, V_t* v2) {              // O(1) operation
        return adjacency_matrix[v1].count(v2) 
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
    add_edge(E_t* e, bool is_undirected=true) {     // O(1) operation
        auto src = (V_t*)e->src;
        auto dst = (V_t*)e->dst;
        if (!contains(src) || !contains(dst))   return false;
        if (contains(src, dst))                 return false;
        edges.push_back(e);
        tlm::put(adjacency_matrix, src, dst, e);
        if (is_undirected)  tlm::put(adjacency_matrix, dst, src, e);
        adjacency_lists[src].push_back(dst);
        if (is_undirected)  adjacency_lists[dst].push_back(src);
        graph_has_changed = true;
        return true;
    }

    virtual V_t*
    get_vertex(uint id) {        // O(1) operation
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

        for (auto w : adjacency_lists[v]) {
            auto& adj = adjacency_lists[w];
            for (auto it = adj.begin(); it != adj.end();) {
                if (*it == v)   it = adj.erase(it);
                else            it++;
            }
            adjacency_matrix[v][w] = nullptr;
            adjacency_matrix[w][v] = nullptr;
        }
        adjacency_lists.erase(v);
        for (auto it = edges.begin(); it != edges.end();) {
            auto e = *it;
            auto u1 = (V_t*)e->src;
            auto u2 = (V_t*)e->dst;
            if (u1 == v || u2 == v) it = edges.erase(it);
            else                    it++;
        }
        graph_has_changed = true;
        delete v;
    }

    virtual void
    delete_edge(E_t* e) {       // O(m) operation
        if (!contains(e))   return;

        auto src = (V_t*)e->src;
        auto dst = (V_t*)e->dst;

        adjacency_matrix[src][dst] = nullptr;
        adjacency_matrix[dst][src] = nullptr;
        
        auto& adj_src = adjacency_lists[src];
        for (auto it = adj_src.begin(); it != adj_src.end();) {
            if (*it == dst)  it = adj_src.erase(it);
            else                it++;
        }
        auto& adj_dst = adjacency_lists[dst];
        for (auto it = adj_dst.begin(); it != adj_dst.end();) {
            if (*it == src) it = adj_dst.erase(it);
            else            it++;
        }

        for (auto it = edges.begin(); it != edges.end();) {
            if (*it == e)   it = edges.erase(it);
            else            it++;
        }

        graph_has_changed=true;
        delete e;
    }

    std::vector<V_t*>   get_vertices(void) { return vertices; }
    std::vector<E_t*>   get_edges(void) { return edges; }
    std::vector<V_t*>   get_neighbors(V_t* v) { return adjacency_lists[v]; }
    uint                get_degree(V_t* v) { return get_neighbors(v).size(); }
protected:
    std::vector<V_t*>   vertices;
    std::vector<E_t*>   edges;

    TwoLevelMap<V_t*, V_t*, E_t*>       adjacency_matrix;
    std::map<V_t*, std::vector<V_t*>>   adjacency_lists;
    std::map<uint, V_t*>                id_to_vertex;

    bool    graph_has_changed;  // Tracks if the graph has changed. May be useful
                                // for subclasses that need to track state.
};

// Evaluation functions:

namespace search {  // For BFS/DFS

template <class V_t>
using callback_t = std::function<void(V_t*, V_t*)>;    // This callback is
                                                            // called on a 
                                                            // successful traversal
                                                            // from a vertex
                                                            // to its unvisited
                                                            // neighbor

}   // search

// A function that can perform either BFS or DFS (code is essentially
// the same, except for the data structure used to store the vertices:
// queue vs stack, both of which can be modelled by a deque).
template <class V_t, class E_t> void
xfs(Graph<V_t, E_t>* graph, V_t* start, search::callback_t<V_t> cb, bool dfs) {
    std::deque<V_t*> dq;
    std::set<V_t*> visited;
    dq.push_back(start);

    while (dq.size()) {
        V_t* v;
        if (dfs) {
            v = dq.back();
            dq.pop_back();
        } else {
            v = dq.front();
            dq.pop_front();
        }

        for (auto w : graph->get_neighbors(v)) {
            if (visited.count(w))   continue;
            cb(v, w);
            dq.push_back(w);
        }
        visited.insert(v);
    }
}

namespace dijkstra {
    template <class V_t>
    using ewf_t = std::function<fp_t(V_t*, V_t*)>;

    template <class V_t, class data_t>
    using DistanceMatrix = TwoLevelMap<V_t*, V_t*, data_t>;

    // The callback below is for creating a distance matrix over multiple
    // calls of dijkstra's. The function should return some data for a
    // matrix entry given (1, 2) two vertices indexing into the matrix,
    // (3) the distances between the first input and every other vertex,
    // (4) the predecessor entry for each vertex in a path from the 
    // first input to that vertex.
    template <class V_t, class data_t>
    using callback_t = std::function<data_t(V_t*, V_t*, 
                            const std::map<V_t*, fp_t>&, const std::map<V_t*, V_t*>&)>;
}   // dijkstra

// Performs Dijkstra's algorithm given an edge weight function (ewf_t).
template <class V_t, class E_t> void
exec_dijkstra(
    Graph<V_t, E_t>* graph, 
    V_t* src,
    std::map<V_t*, fp_t>& distances,
    std::map<V_t*, V_t*>& predecessors,
    dijkstra::ewf_t<V_t> edge_w_func)
{
    typedef struct { V_t* v; fp_t s; } pqv_t;
    struct cmp {
        bool operator()(const pqv_t& v1, const pqv_t& v2) {
            return v1.s > v2.s;
        }
    };

    std::map<V_t*, pqv_t> v2pv;
    std::priority_queue<pqv_t, std::vector<pqv_t>, cmp> queue;
    for (V_t* v : graph->get_vertices()) {
        if (v == src)   distances[v] = 0;
        else            distances[v] = std::numeric_limits<fp_t>::max();
        predecessors[v] = v;

        pqv_t pv = {v, distances[v]};
        queue.push(pv);
        v2pv[v] = pv;
    }

    while (!queue.empty()) {
        pqv_t pv = queue.top();
        auto v = pv.v;
        queue.pop();
        if (pv.s != distances[v])  continue;   // This entry is outdated.

        auto adj = graph->get_neighbors(v);
        for (V_t* w : adj) {
            fp_t new_dist = distances[v] + edge_w_func(v, w);
            if (new_dist < distances[w]) {
                distances[w] = new_dist;
                predecessors[w] = v;
                pqv_t pw = {w, new_dist};
                queue.push(pw);
            }
        }
    }
}

// Creates a distance matrix by calling dijkstra's multiple times.
// Each entry is populated according to the callback (defined above in the
// header).
template <class V_t, class E_t, class data_t> dijkstra::DistanceMatrix<V_t, data_t>
create_distance_matrix(
    Graph<V_t, E_t>* graph, 
    dijkstra::ewf_t<V_t> edge_w_func,
    dijkstra::callback_t<V_t, data_t> cb) 
{
    dijkstra::DistanceMatrix<V_t, data_t> mat;
    auto vertices = graph->get_vertices();
    for (uint i = 0; i < vertices.size(); i++) {
        V_t* src = vertices[i];
        std::map<V_t*, V_t*> pred;
        std::map<V_t*, fp_t> dist;
        exec_dijkstra(graph, src, dist, pred, edge_w_func);
        for (uint j = 0; j < vertices.size(); j++) {
            V_t* dst = vertices[j];
            data_t x = cb(src, dst, dist, pred);
            tlm::put(mat, src, dst, x);
        }
    }
    return mat;
}


}   // graph
}   // qontra

#endif  // GRAPH_h
