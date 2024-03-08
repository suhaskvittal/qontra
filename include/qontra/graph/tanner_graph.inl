/*
 *  author: Suhas Vittal
 *  date    16 February 2024
 * */

namespace qontra {
namespace graph {

inline bool
TannerGraph::add_vertex(sptr<tanner::vertex_t> v) {
    if (!Graph::add_vertex(v))  return false;
    get_vertices_by_type_(v->qubit_type).push_back(v);
    return true;
}

inline bool
TannerGraph::add_edge(sptr<tanner::edge_t> e) {
    sptr<tanner::vertex_t> v = e->get_source<tanner::vertex_t>(),
                           w = e->get_target<tanner::vertex_t>();
    // Make sure the edge preserves the bipartite property.
    bool src_is_parity = v->qubit_type != tanner::vertex_t::type::data;
    bool dst_is_parity = w->qubit_type != tanner::vertex_t::type::data;
    if (src_is_parity == dst_is_parity) return false;
    return Graph::add_edge(e);
}

inline void
TannerGraph::delete_vertex(sptr<tanner::vertex_t> v) {
    std::vector<sptr<tanner::vertex_t>>& vlist = get_vertices_by_type_(v->qubit_type);
    for (auto it = vlist.begin(); it != vlist.end();) {
        if (*it == v)   it = vlist.erase(it);
        else            it++;
    }
    Graph::delete_vertex(v);
}

inline std::vector<sptr<tanner::vertex_t>>
TannerGraph::get_vertices_by_type(tanner::vertex_t::type t) const {
    if (t == tanner::vertex_t::type::data) return data_qubits;
    else if (t == tanner::vertex_t::type::xparity) return xparity_checks;
    else return zparity_checks;
}

inline std::vector<sptr<tanner::vertex_t>>
TannerGraph::get_checks() const {
    std::vector<sptr<tanner::vertex_t>> parity_qubits(xparity_checks);
    for (auto c : zparity_checks)   parity_qubits.push_back(c);
    return parity_qubits;
}

inline std::vector<TannerGraph::obs_t>
TannerGraph::get_obs(bool get_x_obs) const {
    if (get_x_obs)  return x_obs_list;
    else            return z_obs_list;
}

inline std::vector<sptr<tanner::vertex_t>>&
TannerGraph::get_vertices_by_type_(tanner::vertex_t::type t) {
    if (t == tanner::vertex_t::type::data) return data_qubits;
    else if (t == tanner::vertex_t::type::xparity) return xparity_checks;
    else return zparity_checks;
}

inline int
TannerGraph::compute_check_color_map(std::map<sptr<tanner::vertex_t>, int>& color_map) const {
    std::map<sptr<tanner::vertex_t>, int> x_color_map;
    int max_color_z = update_check_color_map(color_map, false);
    int max_color_x = update_check_color_map(x_color_map, true);
    for (auto& [ v, c ] : x_color_map) color_map[v] = c;
    return std::max(max_color_x, max_color_z);
}

template <> inline std::string
print_v(sptr<tanner::vertex_t> v) {
    std::string type_prefix;
    if (v->qubit_type == tanner::vertex_t::type::data) {
        type_prefix = "d";
    } else if (v->qubit_type == tanner::vertex_t::type::xparity) {
        type_prefix = "x";
    } else {
        type_prefix = "z";
    }
    return type_prefix + std::to_string(v->id & tanner::VERTEX_ID_NUMBER_MASK);
}

}   // graph
}   // qontra
