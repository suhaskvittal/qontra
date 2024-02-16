/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#define __INC_VERTEX__  Hypergraph<V, HE>::inc_vertex_t

namespace qontra {
namespace graph {

namespace base {

template <class V> inline std::vector<sptr<V>>
hyperedge_t::operator()() const {
    return std::vector<sptr<V>>(endpoints.begin(), endpoints.end());
}

template <class V> inline sptr<V>
hyperedge_t::operator[](size_t i) const {
    return std::reinterpret_pointer_cast<V>(endpoints.at(i));
}

inline size_t
hyperedge_t::get_order() const {
    return endpoints.size();
}

}   // base

template <class V, class HE> inline void
Hypergraph<V, HE>::change_id(sptr<V> v, uint64_t to) {
    id_to_vertex.erase(v->id);
    v->id = to;
    id_to_vertex[to] = v;
}

template <class V, class HE> inline void
Hypergraph<V, HE>::manual_update_id(sptr<V> v, uint64_t old_id, uint64_t new_id) {
    id_to_vertex.erase(old_id);
    id_to_vertex[new_id] = v;
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::contains(uint64_t id) const {
    return id_to_vertex.count(id);
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::contains(sptr<V> v) const {
    return id_to_vertex.count(v->id);
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::contains(sptr<HE> e) const {
    sptr<__INC_VERTEX__> iv_e = get_incidence_vertex(e);
    return iv_e != nullptr;
}

template <class V, class HE>
template <class CONTAINER>
inline bool
Hypergraph<V, HE>::contains(CONTAINER vlist) const {
    std::vector<__INC_VERTEX__> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist.push_back(get_incidence_vertex(vlist.at(i)));
    }
    std::vector<__INC_VERTEX__> common = incidence_graph.get_common_neighbors(inc_vlist);
    return common.size() == 1 && incidence_graph.get_degree(common[0]) == vlist.size();
}

template <class V, class HE>
template <class ITER>
inline bool
Hypergraph<V, HE>::contains(ITER begin, ITER end) const {
    return contains(std::vector<V>(begin, end));
}

template <class V, class HE> inline sptr<V>
Hypergraph<V, HE>::make_vertex() const {
    return std::make_shared<V>();
}

template <class V, class HE> inline sptr<V>
Hypergraph<V, HE>::make_vertex(uint64_t id) const {
    sptr<V> v = std::make_shared<V>();
    v->id = id;
    return v;
}

template <class V, class HE>
template <class CONTAINER>
inline sptr<HE>
Hypergraph<V, HE>::make_edge(CONTAINER vlist) const {
    sptr<HE> e = std::make_shared<HE>();
    e->endpoints = std::vector<sptr<void>>(vlist.begin(), vlist.end());
    return e;
}

template <class V, class HE>
template <class ITER>
inline sptr<HE>
Hypergraph<V, HE>::make_edge(ITER begin, ITER end) const {
    return make_edge(std::vector<V>(begin, end));
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::add_vertex(sptr<V> v) {
    if (v == nullptr || contains(v->id)) return false;
    id_to_vertex[v->id] = v;
    vertices.push_back(v);
    // Add to incidence graph.
    incidence_graph.make_and_add_vertex(reinterpret_cast<uint64_t>(v));
    graph_has_changed = true;
    return true;
}

template <class V, class HE> bool
Hypergraph<V, HE>::add_edge(sptr<HE> e) {
    if (!is_edge_valid(e) || contains(e.endpoints.begin(), e.endpoints.end())) {
        return false;
    }
    edges.push_back(e);
    // Update incidence_graph and adjacency_lists
    sptr<__INC_VERTEX__> iv_e = incidence_graph.make_and_add_vertex(reinterpret_cast<uint64_t>(e));
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> v = e[i];
        // Update adjacency_list before graph as share_hyperedge uses the
        // incidence_graph.
        for (size_t j = i+1; j < e->get_order(); j++) {
            sptr<V> w = e[j];
            if (!share_hyperedge(v, w)) {
                adjacency_lists[v].push_back(w);
                vtils::tlm_put(adjacency_mult_map, v, w, 0);
            }
            adjacency_mult_map[v][w]++;
        }
        sptr<__INC_VERTEX__> iv_v = get_incidence_vertex(v);
        incidence_graph.make_and_add_edge(iv_v, iv_e);
    }
    graph_has_changed = true;
    return true;
}

template <class V, class HE> inline sptr<V>
Hypergraph<V, HE>::make_and_add_vertex(uint64_t id) {
    sptr<V> v = make_vertex(id);
    return add_vertex(v) ? v : nullptr;
}

template <class V, class HE>
template <class CONTAINER>
inline sptr<HE>
Hypergraph<V, HE>::make_and_add_edge(CONTAINER vlist) {
    sptr<E> e = make_edge(vlist);
    return add_edge(e) ? e : nullptr;
}

template <class V, class HE>
template <class ITER>
inline sptr<HE>
Hypergraph<V, HE>::make_and_add_edge(ITER begin, ITER end) {
    sptr<E> e = make_edge(begin, end);
    return add_edge(e) ? e : nullptr;
}

template <class V, class HE> inline sptr<V>
Hypergraph<V, HE>::get_vertex(uint64_t id) const {
    return id_to_vertex.at(id);
}

template <class V, class HE> 
template <class CONTAINER>
inline sptr<HE>
Hypergraph<V, HE>::get_edge(CONTAINER vlist) const {
    std::vector<__INC_VERTEX__> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist.push_back(get_incidence_vertex(vlist.at(i)));
    }
    std::vector<__INC_VERTEX__> common = incidence_graph.get_common_neighbors(inc_vlist);
    return common.size() == 1 && incidence_graph.get_degree(common[0]) == vlist.size() ? common[0] : nullptr;
}

template <class V, class HE>
template <class ITER>
inline sptr<HE>
Hypergraph<V, HE>::get_edge(ITER begin, ITER end) const {
    return get_edge(std::vector<V>(begin, end));
}

template <class V, class HE>
template <class CONTAINER>
inline sptr<HE>
Hypergraph<V, HE>::get_edge_from_ids(CONTAINER idlist) const {
    std::vector<sptr<V>> vlist(idlist.size());
    for (size_t i = 0; i < idlist.size(); i++) {
        vlist[i] = get_vertex(idlist.at(i));
    }
    return get_edge(vlist);
}

template <class V, class HE>
template <class ITER>
inline sptr<HE>
Hypergraph<V, HE>::get_edge_from_ids(ITER begin, ITER end) const {
    return get_edge_from_ids(std::vector<uint64_t>(begin, end));
}

template <class V, class HE> void
Hypergraph<V, HE>::delete_vertex(sptr<V> v) {
    if (v == nullptr || !contains(v)) return;
    for (auto it = vertices.begin(); it != vertices.end(); ) {
        if (*it == v)   it = vertices.erase(it);
        else            it++;
    }
    adjacency_lists.erase(v);
    // Delete v from the incidence graph.
    sptr<__INC_VERTEX__> iv_v = get_incidence_vertex(v);
    incidence_graph.delete_vertex(iv_v);
    // Delete edges containing v.
    for (auto it = edges.begin(); it != edges.end(); ){
        sptr<HE> e = *it;
        if (has_endpoint(e, v)) {
            // Delete from the incidence_graph.
            sptr<__INC_VERTEX__> iv_e = get_incidence_vertex(e);
            incidence_graph.delete_vertex(iv_e);
            // Update the adjacency lists of all endpoints of e.
            for (size_t i = 0; i < e->get_order(); i++) {
                sptr<V> w = e[i];
                if (v == w) continue;
                auto& adj = adjacency_lists[w];
                for (auto iit = adj.begin(); iit != adj.end(); ) {
                    if (*iit == v) {
                        adjacency_mult_map[w].erase(*iit);
                        it = adj.erase(it);
                    } else {
                        if ((--adjacency_mult_map[w][*iit]) == 0) {
                            adjacency_mult_map[w].erase(*iit);
                            it = adj.erase(it);
                        } else {
                            it++;
                        }
                    }
                }
            }
            it = edges.erase(it);
        } else {
            it++;
        }
    }
    graph_has_changed = true;
}

template <class V, class HE> void
Hypergraph<V, HE>::delete_edge(sptr<HE> e ) {
    if (e == nullptr || !contains(e)) return;
    // Delete from the incidence_graph.
    sptr<__INC_VERTEX__> iv_e = get_incidence_vertex(e);
    incidence_graph.delete_vertex(iv_e);
    // Update the adjacency lists of all endpoints of e.
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> w = e[i];
        if (v == w) continue;
        auto& adj = adjacency_lists[w];
        for (auto iit = adj.begin(); iit != adj.end(); ) {
            if (*iit == v) {
                adjacency_mult_map[w].erase(*iit);
                it = adj.erase(it);
            } else {
                if ((--adjacency_mult_map[w][*iit]) == 0) {
                    adjacency_mult_map[w].erase(*iit);
                    it = adj.erase(it);
                } else {
                    it++;
                }
            }
        }
    }
    for (auto it = edges.begin(); it != edges.end();) {
        if (*it == e)   it = edges.erase(it);
        else            it++;
    }
    graph_has_changed = true;
}

template <class V, class HE> inline std::vector<sptr<V>>
Hypergraph<V, HE>::get_vertices() const {
    return vertices;
}

template <class V, class HE> inline std::vector<sptr<HE>>
Hypergraph<V, HE>::get_edges() const {
    return edges;
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::n() const {
    return vertices.size();
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::m() const {
    return edges.size();
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::has_endpoint(sptr<HE> e, sptr<V> v) const {
    sptr<__INC_VERTEX__> iv_e = get_incidence_vertex(e),
                         iv_v = get_incidence_vertex(v);
    return incidence_graph.contains(iv_e, iv_v);
}

template <class V, class HE>
template <class CONTAINER>
inline bool
Hypergraph<V, HE>::share_hyperedge(CONTAINER vlist) const {
    std::vector<sptr<__INC_VERTEX__>> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist[i] = get_incidence_vertex(vlist.at(i));
    }
    return incidence_graph.get_common_neighbors(inc_vlist).size();
}

template <class V, class HE> inline std::vector<sptr<V>>
Hypergraph<V, HE>::get_neighbors(sptr<V> v) const {
    return adjacency_lists.at(v);
}

template <class V, class HE>
template <class CONTAINER>
inline std::vector<sptr<V>>
Hypergraph<V, HE>::get_common_neighbors(CONTAINER vlist) const {
    if (vlist.empty()) return {};
    sptr<V> common = get_neighbors(vlist.at(0));
    for (size_t i = 1; i < vlist.size(); i++) {
        for (auto it = common.begin(); it != common.end(); ) {
            if (!share_hyperedge(vlist.at(i), *it)) it = common.erase(it);
            else                                    it++;
        }
    }
    return common;
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::get_degree(sptr<V> v) const {
    sptr<__INC_VERTEX__> iv_v = get_incidence_vertex(v);
    return incidence_graph.get_degree(iv_v);
}

template <class V, class HE> inline fp_t
Hypergraph<V, HE>::get_mean_degree() {
    update_state();
    return mean_degree;
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::get_max_degree() {
    update_state();
    return max_degree;
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::get_max_order() {
    update_state();
    return max_order;
}

template <class V, class HE> inline fp_t
Hypergraph<V, HE>::const_get_mean_degree() const {
    if (!graph_has_changed) return mean_degree;
    fp_t sdg = 0.0;
    for (sptr<V> v : vertices) {
        sdg += static_cast<fp_t>(get_degree(v));
    }
    return sdg / static_cast<fp_t>(n());
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::const_get_max_degree() const {
    if (!graph_has_changed) return max_degree;
    size_t mdg = 0;
    for (sptr<V> v : vertices) {
        mdg = std::max(get_degree(v), mdg);
    }
    return mdg;
}

template <class V, class HE> inline size_t
Hypergraph<V, HE>::const_get_max_order() const {
    if (!graph_has_changed) return max_order;
    size_t mord = 0;
    for (sptr<V> e : edges) {
        mord = std::max(e->get_order(), mord);
    }
    return mord;
}

template <class V, class HE> inline void
Hypergraph<V, HE>::force_update_state() {
    graph_has_changed = true;
    update_state();
}

template <class V, class HE> inline bool
Hypergraph<V, HE>::update_state() {
    if (!graph_has_changed) return false;
    mean_degree = const_get_mean_degree();
    max_degree = const_get_max_degree();
    max_order = const_get_max_order();
    graph_has_changed = false;
    return true;
}

template <class V, class HE>
template <class PTR>
inline sptr<__INC_VERTEX__>
Hypergraph<V, HE>::get_incidence_vertex(PTR obj_p) const {
    return incidence_graph.get_vertex(reinterpret_cast<uint64_t>(obj_p));
}

template <class V, class HE> bool
Hypergraph<V, HE>::is_edge_valid(sptr<HE> e) {
    if (e == nullptr) return false;
    // Check that all endpoints are good. That is:
    //  1. None are nullptr.
    //  2. None are repeated.
    //  3. All of them exist in this graph.
    std::set<sptr<V>> visited;
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> v = e[i];
        if (v == nullptr) return false;
        if (visited.count(v)) return false;
        if (!contains(v)) return false;
        visited.insert(v);
    }
    return true;
}

}   // graph
}   // qontra
