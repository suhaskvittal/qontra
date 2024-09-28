/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "qontra/graph/tanner_graph.h"
#include "qontra/graph/algorithms/coloring.h"

#include <limits>

namespace qontra {
namespace graph {

using namespace tanner;

bool
TannerGraph::add_vertex(sptr<vertex_t> v) {
    if (!Graph::add_vertex(v))  return false;
    get_vertices_by_type_(v->qubit_type).push_back(v);
    return true;
}

bool
TannerGraph::add_edge(sptr<edge_t> e) {
    sptr<vertex_t> v = e->get_source<vertex_t>(),
                           w = e->get_target<vertex_t>();
    // Make sure the edge preserves the bipartite property.
    bool src_is_parity = v->qubit_type != vertex_t::type::data;
    bool dst_is_parity = w->qubit_type != vertex_t::type::data;
    if (src_is_parity == dst_is_parity) return false;
    return Graph::add_edge(e);
}

void
TannerGraph::delete_vertex(sptr<vertex_t> v) {
    std::vector<sptr<vertex_t>>& vlist = get_vertices_by_type_(v->qubit_type);
    for (auto it = vlist.begin(); it != vlist.end();) {
        if (*it == v)   it = vlist.erase(it);
        else            it++;
    }
    Graph::delete_vertex(v);
}

std::vector<sptr<vertex_t>>
TannerGraph::get_vertices_by_type(vertex_t::type t) const {
    if (t == vertex_t::type::data) return data_qubits;
    else if (t == vertex_t::type::xparity) return xparity_checks;
    else return zparity_checks;
}

std::vector<sptr<vertex_t>>
TannerGraph::get_checks() const {
    std::vector<sptr<vertex_t>> parity_qubits(xparity_checks);
    for (auto c : zparity_checks)   parity_qubits.push_back(c);
    return parity_qubits;
}

std::vector<TannerGraph::obs_t>
TannerGraph::get_obs(bool get_x_obs) const {
    if (get_x_obs)  return x_obs_list;
    else            return z_obs_list;
}

std::vector<TannerGraph::obs_t>&
TannerGraph::get_obs_ref(bool get_x_obs) {
    if (get_x_obs)  return x_obs_list;
    else            return z_obs_list;
}

std::vector<sptr<vertex_t>>&
TannerGraph::get_vertices_by_type_(vertex_t::type t) {
    if (t == vertex_t::type::data) return data_qubits;
    else if (t == vertex_t::type::xparity) return xparity_checks;
    else return zparity_checks;
}

int
TannerGraph::compute_check_color_map(std::unordered_map<sptr<vertex_t>, int>& color_map) const {
    std::unordered_map<sptr<vertex_t>, int> x_color_map;
    int max_color_z = update_check_color_map(color_map, false);
    int max_color_x = update_check_color_map(x_color_map, true);
    for (auto& [ v, c ] : x_color_map) color_map[v] = c;
    return std::max(max_color_x, max_color_z);
}

int
TannerGraph::update_check_color_map(std::unordered_map<sptr<vertex_t>, int>& check_color_map, bool use_x_checks) const {
    const std::vector<sptr<vertex_t>>& checks = use_x_checks ? xparity_checks : zparity_checks;
    // Create an interaction graph for these checks.
    uptr<Graph<vertex_t, base::edge_t>> gr = std::make_unique<Graph<vertex_t, base::edge_t>>();
    for (sptr<vertex_t> v : checks) gr->add_vertex(v);
    // Add edges if two vertices share neighbors in the Tanner graph.
    for (size_t i = 0; i < checks.size(); i++) {
        sptr<vertex_t> v = checks.at(i);
        for (size_t j = i+1; j < checks.size(); j++) {
            sptr<vertex_t> w = checks.at(j);
            if (get_common_neighbors({v, w}).size()) {
                gr->make_and_add_edge(v, w);
            }
        }
    }
    int best_coloring = k_coloring_rlf(gr.get(), check_color_map);
    for (size_t s = 0; s < checks.size() && (best_coloring > 2 || best_coloring < 0); s++) {
        std::unordered_map<sptr<vertex_t>, int> cm;
        int c = k_coloring_greedy(gr.get(), cm, s);
        if (c < best_coloring) {
            check_color_map = std::move(cm);
            cm.clear();
            best_coloring = c;
        }
    }
    return best_coloring;
}

template <> std::string
print_v(sptr<vertex_t> v) {
    std::string type_prefix;
    if (v->qubit_type == vertex_t::type::data) {
        type_prefix = "d";
    } else if (v->qubit_type == vertex_t::type::xparity) {
        type_prefix = "x";
    } else {
        type_prefix = "z";
    }
    return type_prefix + std::to_string(v->id & VERTEX_ID_NUMBER_MASK);
}

}   // graph
}   // qontra
