/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

namespace qontra {
namespace graph {

namespace base {

template <class V> inline std::vector<sptr<V>>
hyperedge_t::get() const {
    return std::vector<sptr<V>>(endpoints.begin(), endpoints.end());
}

template <class V> inline sptr<V>
hyperedge_t::get(size_t i) const {
    return std::reinterpret_pointer_cast<V>(endpoints.at(i));
}

inline size_t
hyperedge_t::get_order() const {
    return endpoints.size();
}

}   // base

template <class V, class HE> inline std::string
print_he(sptr<HE> e) {
    std::string out = "(";
    for (size_t i = 0; i < e.get_order(); i++) {
        if (i > 0) out += ",";
        out += print_v<V>(e->endpoints[i]);
    }
    out += ")";
    return out;
}

template <class V, class HE> inline void
HyperGraph<V, HE>::change_id(sptr<V> v, uint64_t to) {
    id_to_vertex.erase(v->id);
    v->id = to;
    id_to_vertex[to] = v;
}

template <class V, class HE> inline void
HyperGraph<V, HE>::manual_update_id(sptr<V> v, uint64_t old_id, uint64_t new_id) {
    id_to_vertex.erase(old_id);
    id_to_vertex[new_id] = v;
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::contains(uint64_t id) const {
    return id_to_vertex.count(id);
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::contains(sptr<V> v) const {
    return id_to_vertex.count(v->id);
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::contains(sptr<HE> e) const {
    sptr<base::vertex_t> iv_e = get_incidence_vertex(e);
    return iv_e != nullptr;
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::contains(std::vector<sptr<void>> vlist) const {
    std::vector<sptr<V>> _vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        _vlist[i] = std::reinterpret_pointer_cast<V>(vlist[i]);
    }
    return contains(_vlist);
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::contains(std::vector<sptr<V>> vlist) const {
    std::vector<sptr<base::vertex_t>> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist.push_back(get_incidence_vertex(vlist.at(i)));
    }
    std::vector<sptr<base::vertex_t>> common = incidence_graph.get_common_neighbors(inc_vlist);
    return common.size() == 1 && incidence_graph.get_degree(common[0]) == vlist.size();
}

template <class V, class HE> inline sptr<V>
HyperGraph<V, HE>::make_vertex() const {
    return std::make_shared<V>();
}

template <class V, class HE> inline sptr<V>
HyperGraph<V, HE>::make_vertex(uint64_t id) const {
    sptr<V> v = std::make_shared<V>();
    v->id = id;
    return v;
}

template <class V, class HE> inline sptr<HE>
HyperGraph<V, HE>::make_edge(std::vector<sptr<V>> vlist) const {
    sptr<HE> e = std::make_shared<HE>();
    e->endpoints = std::vector<sptr<void>>(vlist.begin(), vlist.end());
    return e;
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::add_vertex(sptr<V> v) {
    if (v == nullptr || contains(v->id)) return false;
    id_to_vertex[v->id] = v;
    vertices.push_back(v);
    // Add to incidence graph.
    sptr<base::vertex_t> iv_v =
        incidence_graph.make_and_add_vertex(reinterpret_cast<uint64_t>(v.get()));
    incidence_object_map[iv_v] = std::reinterpret_pointer_cast<void>(v);
    graph_has_changed = true;
    return true;
}

template <class V, class HE> bool
HyperGraph<V, HE>::add_edge(sptr<HE> e) {
    if (!is_edge_valid(e) || contains(e->endpoints)) {
        return false;
    }
    edges.push_back(e);
    // Update incidence_graph and adjacency_lists
    sptr<base::vertex_t> iv_e = incidence_graph.make_and_add_vertex(reinterpret_cast<uint64_t>(e.get()));
    incidence_object_map[iv_e] = std::reinterpret_pointer_cast<void>(e);
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> v = std::reinterpret_pointer_cast<V>(e->endpoints[i]);
        // Update adjacency_list before graph as share_hyperedge uses the
        // incidence_graph.
        for (size_t j = i+1; j < e->get_order(); j++) {
            sptr<V> w = std::reinterpret_pointer_cast<V>(e->endpoints[j]);
            if (!share_hyperedge({v, w})) {
                adjacency_lists[v].push_back(w);
                adjacency_lists[w].push_back(v);
                vtils::tlm_put(adjacency_mult_map, v, w, static_cast<size_t>(0));
                vtils::tlm_put(adjacency_mult_map, w, v, static_cast<size_t>(0));
            }
            adjacency_mult_map[v][w]++;
            adjacency_mult_map[w][v]++;
        }
        sptr<base::vertex_t> iv_v = get_incidence_vertex(v);
        incidence_graph.make_and_add_edge(iv_v, iv_e);
    }
    graph_has_changed = true;
    return true;
}

template <class V, class HE> inline sptr<V>
HyperGraph<V, HE>::make_and_add_vertex(uint64_t id) {
    sptr<V> v = make_vertex(id);
    return add_vertex(v) ? v : nullptr;
}

template <class V, class HE> inline sptr<HE>
HyperGraph<V, HE>::make_and_add_edge(std::vector<sptr<V>> vlist) {
    sptr<HE> e = make_edge(vlist);
    return add_edge(e) ? e : nullptr;
}

template <class V, class HE> inline sptr<V>
HyperGraph<V, HE>::get_vertex(uint64_t id) const {
    return id_to_vertex.at(id);
}

template <class V, class HE>  inline sptr<HE>
HyperGraph<V, HE>::get_edge(std::vector<sptr<V>> vlist) const {
    std::vector<sptr<base::vertex_t>> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist.push_back(get_incidence_vertex(vlist.at(i)));
    }
    std::vector<sptr<base::vertex_t>> common = incidence_graph.get_common_neighbors(inc_vlist);
    return common.size() == 1 && incidence_graph.get_degree(common[0]) == vlist.size()
            ? std::reinterpret_pointer_cast<HE>(incidence_object_map.at(common[0])) : nullptr;
}

template <class V, class HE> inline sptr<HE>
HyperGraph<V, HE>::get_edge(std::vector<uint64_t> idlist) const {
    std::vector<sptr<V>> vlist(idlist.size());
    for (size_t i = 0; i < idlist.size(); i++) {
        vlist[i] = get_vertex(idlist.at(i));
    }
    return get_edge(vlist);
}

template <class V, class HE> void
HyperGraph<V, HE>::delete_vertex(sptr<V> v) {
    if (v == nullptr || !contains(v)) return;
    for (auto it = vertices.begin(); it != vertices.end(); ) {
        if (*it == v)   it = vertices.erase(it);
        else            it++;
    }
    for (sptr<V> w : adjacency_lists[v]) {
        auto& adj = adjacency_lists[w];
        for (auto it = adj.begin(); it != adj.end(); ) {
            if (*it == v)   it = adj.erase(it);
            else            it++;
        }
        adjacency_mult_map[w].erase(v);
    }
    adjacency_lists.erase(v);
    // Delete v from the incidence graph.
    sptr<base::vertex_t> iv_v = get_incidence_vertex(v);
    incidence_graph.delete_vertex(iv_v);
    // Delete edges containing v.
    for (auto it = edges.begin(); it != edges.end(); ){
        sptr<HE> e = *it;
        if (has_endpoint(e, v)) {
            // Delete from the incidence_graph.
            sptr<base::vertex_t> iv_e = get_incidence_vertex(e);
            incidence_graph.delete_vertex(iv_e);
            update_adjacency_lists_after_delete(e);
            it = edges.erase(it);
        } else {
            it++;
        }
    }
    graph_has_changed = true;
}

template <class V, class HE> void
HyperGraph<V, HE>::delete_edge(sptr<HE> e) {
    if (e == nullptr || !contains(e)) return;
    // Delete from the incidence_graph.
    sptr<base::vertex_t> iv_e = get_incidence_vertex(e);
    incidence_graph.delete_vertex(iv_e);
    // Update the adjacency lists of all endpoints of e.
    update_adjacency_lists_after_delete(e);
    for (auto it = edges.begin(); it != edges.end();) {
        if (*it == e)   it = edges.erase(it);
        else            it++;
    }
    graph_has_changed = true;
}

template <class V, class HE> inline std::vector<sptr<V>>
HyperGraph<V, HE>::get_vertices() const {
    return vertices;
}

template <class V, class HE> inline std::vector<sptr<HE>>
HyperGraph<V, HE>::get_edges() const {
    return edges;
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::n() const {
    return vertices.size();
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::m() const {
    return edges.size();
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::has_endpoint(sptr<HE> e, sptr<V> v) const {
    sptr<base::vertex_t> iv_e = get_incidence_vertex(e),
                         iv_v = get_incidence_vertex(v);
    return incidence_graph.contains(iv_e, iv_v);
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::share_hyperedge(std::vector<sptr<V>> vlist) const {
    std::vector<sptr<base::vertex_t>> inc_vlist(vlist.size());
    for (size_t i = 0; i < vlist.size(); i++) {
        inc_vlist[i] = get_incidence_vertex(vlist.at(i));
    }
    return incidence_graph.get_common_neighbors(inc_vlist).size();
}

template <class V, class HE> inline std::vector<sptr<V>>
HyperGraph<V, HE>::get_neighbors(sptr<V> v) const {
    return adjacency_lists.at(v);
}

template <class V, class HE> inline std::vector<sptr<V>>
HyperGraph<V, HE>::get_common_neighbors(std::vector<sptr<V>> vlist) const {
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
HyperGraph<V, HE>::get_degree(sptr<V> v) const {
    sptr<base::vertex_t> iv_v = get_incidence_vertex(v);
    return incidence_graph.get_degree(iv_v);
}

template <class V, class HE> inline fp_t
HyperGraph<V, HE>::get_mean_degree() {
    update_state();
    return mean_degree;
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::get_max_degree() {
    update_state();
    return max_degree;
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::get_max_order() {
    update_state();
    return max_order;
}

template <class V, class HE> inline fp_t
HyperGraph<V, HE>::const_get_mean_degree() const {
    if (!graph_has_changed) return mean_degree;
    fp_t sdg = 0.0;
    for (sptr<V> v : vertices) {
        sdg += static_cast<fp_t>(get_degree(v));
    }
    return sdg / static_cast<fp_t>(n());
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::const_get_max_degree() const {
    if (!graph_has_changed) return max_degree;
    size_t mdg = 0;
    for (sptr<V> v : vertices) {
        mdg = std::max(get_degree(v), mdg);
    }
    return mdg;
}

template <class V, class HE> inline size_t
HyperGraph<V, HE>::const_get_max_order() const {
    if (!graph_has_changed) return max_order;
    size_t mord = 0;
    for (sptr<HE> e : edges) {
        mord = std::max(e->get_order(), mord);
    }
    return mord;
}

template <class V, class HE> inline void
HyperGraph<V, HE>::force_update_state() {
    graph_has_changed = true;
    update_state();
}

template <class V, class HE> inline bool
HyperGraph<V, HE>::update_state() {
    if (!graph_has_changed) return false;
    mean_degree = const_get_mean_degree();
    max_degree = const_get_max_degree();
    max_order = const_get_max_order();
    graph_has_changed = false;
    return true;
}

template <class V, class HE>
template <class PTR>
inline sptr<base::vertex_t>
HyperGraph<V, HE>::get_incidence_vertex(PTR obj_p) const {
    return incidence_graph.get_vertex(reinterpret_cast<uint64_t>(obj_p.get()));
}

template <class V, class HE> bool
HyperGraph<V, HE>::is_edge_valid(sptr<HE> e) {
    if (e == nullptr) return false;
    // Check that all endpoints are good. That is:
    //  1. None are nullptr.
    //  2. None are repeated.
    //  3. All of them exist in this graph.
    std::set<sptr<V>> visited;
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> v = std::reinterpret_pointer_cast<V>(e->endpoints[i]);
        if (v == nullptr) return false;
        if (visited.count(v)) return false;
        if (!contains(v)) return false;
        visited.insert(v);
    }
    return true;
}

template <class V, class HE> void
HyperGraph<V, HE>::update_adjacency_lists_after_delete(sptr<HE> e) {
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<V> v = std::reinterpret_pointer_cast<V>(e->endpoints[i]);
        auto& adj = adjacency_lists[v];
        for (size_t j = 0; j < e->get_order(); j++) {
            if (i == j) continue;
            sptr<V> w = std::reinterpret_pointer_cast<V>(e->endpoints[j]);
            for (auto it = adj.begin(); it != adj.end(); ) {
                if (*it == w) {
                    if ((--adjacency_mult_map[v][w]) == 0) {
                        adjacency_mult_map[v].erase(w);
                        it = adj.erase(it);
                    } else {
                        it++;
                    }
                }
            }
        }
    }
}

}   // graph
}   // qontra
