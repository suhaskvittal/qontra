/*
 *  author: Suhas Vittal
 *  date:   22 January 2024
 * */

#include <algorithm>

namespace qontra {
namespace graph {

namespace base {

template <class V> inline sptr<V>
edge_t::get_source() {
    return std::reinterpret_pointer_cast<V>(src);
}

template <class V> inline sptr<V>
edge_t::get_target() {
    return std::reinterpret_pointer_cast<V>(dst);
}

}   // base

template <class V> inline std::string
print_v(sptr<V> v) {
    return std::to_string(v->id);
}

template <> inline std::string
print_v(sptr<void> v) {
    return print_v(std::reinterpret_pointer_cast<base::vertex_t>(v));
}

template <class V, class E> inline std::string
print_e(sptr<E> e) {
    sptr<V> src = std::reinterpret_pointer_cast<V>(e->src),
            dst = std::reinterpret_pointer_cast<V>(e->dst);
    std::string dir_prefix = e->is_undirected ? "u" : "d";
    return dir_prefix + "(" + print_v<V>(src) + "," + print_v<V>(dst) + ")";
}


template <class V, class E> inline void
Graph<V, E>::change_id(sptr<V> v, uint64_t to) {
    id_to_vertex.erase(v->id);
    v->id = to;
    id_to_vertex[to] = v;
}

template <class V, class E> inline void
Graph<V, E>::manual_update_id(sptr<V> v, uint64_t old_id, uint64_t new_id) {
    // Use change_id for most cases. This function is only useful
    // if the same vertex is shared across multiple graphs.
    id_to_vertex.erase(old_id);
    id_to_vertex[new_id] = v;
}

template <class V, class E> inline bool
Graph<V, E>::contains(uint64_t id) const {
    return id_to_vertex.count(id);
}

template <class V, class E> inline bool
Graph<V, E>::contains(sptr<V> v) const {
    return id_to_vertex.count(v->id);
}

template <class V, class E> inline bool
Graph<V, E>::contains(sptr<E> e) const {
    sptr<V> v1 = std::reinterpret_pointer_cast<V>(e->src),
            v2 = std::reinterpret_pointer_cast<V>(e->dst);
    return adjacency_matrix.count(v1)
            && adjacency_matrix.at(v1).count(v2)
            && adjacency_matrix.at(v1).at(v2) == e;
}

template <class V, class E> inline bool
Graph<V, E>::contains(sptr<V> v1, sptr<V> v2) const {
    return adjacency_matrix.count(v1)
            && adjacency_matrix.at(v1).count(v2) 
            && (adjacency_matrix.at(v1).at(v2) != nullptr);
}

template <class V, class E> inline sptr<V>
Graph<V, E>::make_vertex() const {
    return std::make_shared<V>();
}

template <class V, class E> inline sptr<V>
Graph<V, E>::make_vertex(uint64_t id) const {
    sptr<V> v = std::make_shared<V>();
    v->id = id;
    return v;
}

template <class V, class E> inline sptr<E>
Graph<V, E>::make_edge(sptr<V> src, sptr<V> dst, bool is_undirected) const {
    sptr<E> e = std::make_shared<E>();
    e->src = src;
    e->dst = dst;
    e->is_undirected = is_undirected;
    return e;
}

template <class V, class E> inline bool
Graph<V, E>::add_vertex(sptr<V> v) {
    if (v == nullptr || contains(v->id)) return false;

    id_to_vertex[v->id] = v;
    vertex_enum_map[v] = vertices.size();
    vertices.push_back(v);
    adjacency_lists[v] = {};

    graph_has_changed = true;
    return true;
}

template <class V, class E> bool
Graph<V, E>::add_edge(sptr<E> e) {
    if (e == nullptr) return false;
    sptr<V> src = std::reinterpret_pointer_cast<V>(e->src),
            dst = std::reinterpret_pointer_cast<V>(e->dst);
    if (src == dst)                         return false;
    if (!contains(src) || !contains(dst))   return false;
    if (contains(src, dst))                 return false;
    edges.push_back(e);

    adjacency_matrix[src][dst] = e;
    if (e->is_undirected) {
        adjacency_matrix[dst][src] = e;
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

template <class V, class E> inline sptr<V>
Graph<V, E>::make_and_add_vertex(uint64_t id) {
    sptr<V> v = make_vertex(id);
    return add_vertex(v) ? v : nullptr;
}

template <class V, class E> inline sptr<E>
Graph<V, E>::make_and_add_edge(sptr<V> v, sptr<V> w, bool is_undirected) {
    sptr<E> e = make_edge(v, w, is_undirected);
    return add_edge(e) ? e : nullptr;
}

template <class V, class E> inline sptr<V>
Graph<V, E>::get_vertex(uint64_t id) const {        // O(1) operation
    if (!id_to_vertex.count(id)) return nullptr;
    return id_to_vertex.at(id);
}

template <class V, class E> inline sptr<E>
Graph<V, E>::get_edge(sptr<V> v1, sptr<V> v2) const {    // O(1) operation
    if (!adjacency_matrix.count(v1) || !adjacency_matrix.at(v1).count(v2)) return nullptr;
    return adjacency_matrix.at(v1).at(v2);
}

template <class V, class E> inline sptr<E>
Graph<V, E>::get_edge(uint64_t id1, uint64_t id2) const {
    return get_edge(get_vertex(id1), get_vertex(id2));
}

template <class V, class E> void
Graph<V, E>::delete_vertex(sptr<V> v) {         // O(n) operation
    if (v == nullptr || !contains(v)) return;
    bool found = false;
    for (auto it = vertices.begin(); it != vertices.end();) {
        if (*it == v) {
            it = vertices.erase(it);
            found = true;
        } else {
            vertex_enum_map[*it] -= static_cast<int>(found);
            it++;
        }
    }
    vertex_enum_map.erase(v);

    for (auto it = edges.begin(); it != edges.end();) {
        sptr<E> e = *it;
        sptr<V> u1 = std::reinterpret_pointer_cast<V>(e->src),
                u2 = std::reinterpret_pointer_cast<V>(e->dst);
        if (u1 == v || u2 == v) { 
            // Delete v from the adjacency list of the other vertex.
            sptr<V> other = u1;
            if (u1 == v)    other = u2;

            auto& adj = adjacency_lists[other];
            for (auto iit = adj.begin(); iit != adj.end();) {
                if (*iit == v)   iit = adj.erase(iit);
                else            iit++;
            }

            auto& r_adj = r_adjacency_lists[other];
            for (auto iit = r_adj.begin(); iit != r_adj.end();) {
                if (*iit == v)   iit = r_adj.erase(iit);
                else            iit++;
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

template <class V, class E> void
Graph<V, E>::delete_edge(sptr<E> e) {  // O(m) operation
    if (e == nullptr || !contains(e)) return;
    sptr<V> src = std::reinterpret_pointer_cast<V>(e->src),
            dst = std::reinterpret_pointer_cast<V>(e->dst);

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

template <class V, class E> inline std::vector<sptr<V>>
Graph<V, E>::get_vertices() const {
    return vertices;
}

template <class V, class E> inline std::vector<sptr<E>>
Graph<V, E>::get_edges() const {
    return edges;
}

template <class V, class E> inline const std::map<sptr<V>, size_t>&
Graph<V, E>::get_enumeration_map() const {
    return vertex_enum_map;
}

template <class V, class E> inline size_t
Graph<V, E>::n() const { 
    return vertices.size();
}

template <class V, class E> inline size_t
Graph<V, E>::m() const { 
    return edges.size();
}

template <class V, class E> inline std::vector<sptr<V>>
Graph<V, E>::get_neighbors(sptr<V> v) const {
    return adjacency_lists.at(v); 
}

template <class V, class E> inline std::vector<sptr<V>>
Graph<V, E>::get_incoming(sptr<V> v) const {
    return r_adjacency_lists.at(v);
}

template <class V, class E> inline std::vector<sptr<V>>
Graph<V, E>::get_outgoing(sptr<V> v) const {
    return adjacency_lists.at(v);
}

template <class V, class E> inline std::vector<sptr<V>>
Graph<V, E>::get_common_neighbors(std::vector<sptr<V>> vlist) const {
    if (vlist.empty()) return {};
    std::vector<sptr<V>> common = get_neighbors(vlist[0]);

    for (size_t i = 1; i < vlist.size(); i++) {
        sptr<V> v = vlist.at(i);
        for (auto it = common.begin(); it != common.end(); ) {
            if (!contains(*it, v))  it = common.erase(it);
            else                    it++;
        }
    }
    return common;
}

template <class V, class E> inline size_t
Graph<V, E>::get_degree(sptr<V> v) const {
    return get_neighbors(v).size();
}

template <class V, class E> inline size_t
Graph<V, E>::get_indegree(sptr<V> v) const {
    return get_incoming(v).size();
}

template <class V, class E> inline size_t
Graph<V, E>::get_outdegree(sptr<V> v) const {
    return get_degree(v);
}

template <class V, class E> inline size_t
Graph<V, E>::get_inoutdegree(sptr<V> v) const {
    return get_indegree(v) + get_outdegree(v);
}

template <class V, class E> inline fp_t
Graph<V, E>::get_mean_degree() const {
    return 2 * static_cast<fp_t>(m()) / static_cast<fp_t>(n()); 
}

template <class V, class E> inline size_t
Graph<V, E>::get_max_degree() {
    update_state();
    return max_degree;
}

template <class V, class E> inline size_t
Graph<V, E>::const_get_max_degree() const {
    if (!graph_has_changed) return max_degree;
    size_t mdg = 0;
    for (sptr<V> v : vertices) mdg = std::max(mdg, get_degree(v));
    return mdg;
}

template <class V, class E> inline void
Graph<V, E>::force_update_state(void) {
    graph_has_changed = true; 
    update_state(); 
}

template <class V, class E> inline bool
Graph<V, E>::update_state() {
    if (!graph_has_changed) return false;
    max_degree = const_get_max_degree();
    graph_has_changed = false;
    return true;
}

}   // graph
}   // qontra
