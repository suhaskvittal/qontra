/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include <algorithm>

#include <lemon/planarity.h>

namespace qontra {
namespace protean {

namespace net {

inline void
phys_vertex_t::consume(sptr<phys_vertex_t> other) {
    for (auto& pair : other->cycle_role_map) {
        size_t cyc = pair.first;
        sptr<raw_vertex_t> rx = pair.second;
        if (cycle_role_map.count(rx)) continue;
        push_back_role(rx, cyc);
    }
}

inline void
phys_vertex_t::add_role(sptr<raw_vertex_t> r, size_t cycle) {
    cycle_role_map.put(r, cycle);
    role_type_set.insert(r->qubit_type);
}

inline void
phys_vertex_t::push_back_role(sptr<raw_vertex_t> r, size_t min_cycle) {
    // Find next available cycle.
    size_t cycle = min_cycle;
    while (true) {
        if (!cycle_role_map.count(cycle)) break;
        cycle++;
    }
    add_role(r, cycle);
}

inline bool
phys_vertex_t::has_role_of_type(raw_vertex_t::type t) {
    return role_type_set.count(t);
}

}   // net

inline sptr<net::raw_vertex_t>
RawNetwork::add_proxy(sptr<net::raw_vertex_t> q1, sptr<net::raw_vertex_t> q2) {
    sptr<raw_edge_t> e = get_edge(q1, q2);
    return add_proxy(e);
}

inline sptr<net::raw_vertex_t>
RawNetwork::proxy_walk(sptr<net::raw_vertex_t> from,
                        sptr<net::raw_vertex_t> thru,
                        std::vector<sptr<net::raw_vertex_t>>& walk_res_ref) 
{
    sptr<net::raw_vertex_t> prev = from, curr = thru;
    while (curr->qubit_type == net::raw_vertex_t::type::proxy) {
        walk_res_ref.push_back(curr);

        sptr<net::raw_vertex_t> next = proxy_indirection_map[thru][curr];
        prev = curr;
        curr = next;
    }
    return curr;
}

inline bool
PhysicalNetwork::add_vertex(sptr<net::phys_vertex_t> v) {
    if (graph::Graph::add_vertex(v)) {
        auto lemon_node = planar_repr.addNode();
        v_phys_lemon_map.put(v, lemon_node);
        return true;
    } else {
        return false;
    }
}

inline bool
PhysicalNetwork::add_edge(sptr<net::phys_edge_t> e) {
    if (graph::Graph::add_edge(e)) {
        if (!test_and_move_edge_to_bulk(e, true)) {
            e->is_out_of_plane = true;
            const size_t k = determine_edge_tsv_layer(e);
            e->tsv_layer = k;
            occupied_tsvs[src].insert(k);
            occupied_tsvs[dst].insert(k);
        }
        return true;
    } else {
        return false;
    }
}

inline void
PhysicalNetwork::delete_vertex(sptr<net::phys_vertex_t> v) {
    // First, erase all traces of v in the tracking structures.
    auto lemon_node = v_phys_lemon_map.at(v);
    v_phys_lemon_map.erase(v);
    for (sptr<net::phys_edge_t> e : get_neighbors(v)) {
        if (e->is_out_of_plane) continue;
        auto lemon_edge = e_phys_lemon_map.at(e);
        e_phys_lemon_map.erase(e);
        planar_repr.erase(lemon_edge);
    }
    // Then delete lemon_node and v.
    planar_repr.erase(lemon_node);
    graph::Graph::delete_vertex(v);
}

inline void
PhysicalNetwork::delete_edge(sptr<net::phys_edge_t> e) {
    if (!e->is_out_of_plane) {
        auto lemon_edge = e_phys_lemon_map.at(e);
        e_phys_lemon_map.erase(e);
        planar_repr.erase(lemon_edge);
    }
    graph::Graph::delete_edge(e);
}

inline bool
PhysicalNetwork::test_and_move_edge_to_bulk(sptr<net::phys_edge_t> e, bool is_new_edge) {
    sptr<net::phys_vertex_t> src = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->src),
                            dst = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->dst);
    auto lemon_src = v_phys_lemon_map.at(src),
            lemon_dst = v_phys_lemon_map.at(dst);
    auto lemon_edge = planar_repr.addEdge(lemon_src, lemon_dst);

    if (is_planar()) {
        e_phys_lemon_map.put(e, lemon_edge);

        if (!is_new_edge) {
            occupied_tsvs[src].erase(e->tsv_layer);
            occupied_tsvs[dst].erase(e->tsv_layer);
        }

        e->is_out_of_plane = false;
        return true;
    } else {
        planar_repr.erase(lemon_edge);
        return false;
    }
}

inline bool
PhysicalNetwork::update_state() {
    if (!graph_has_changed) return false;

    // Update bulk_degree_map.
    for (sptr<net::phys_vertex_t> v : get_vertices()) {
        size_t dg = get_degree(v);
        for (sptr<net::phys_vertex_t> w : get_neighbors(v)) {
            sptr<net::phys_edge_t> e = get_edge(v, w);
            if (e->is_out_of_plane) --dg;
        }
        bulk_degree_map[v] = dg;
    }
}

inline void
PhysicalNetwork::reallocate_edges() {
    // Get all out-of-plane edges, and sort them in order of the maximum degrees of their endpoints.
    // The idea is that if the endpoints have low degree, then an edge is less likely to be out-of-plane.
    std::vector<sptr<net::phys_edge_t>> out_of_plane_edges;
    for (sptr<net::phys_edge_t> e : get_edges()) {
        if (e->is_out_of_plane) out_of_plane_edges.push_back(e);
    }

    std::sort(out_of_plane_edges.begin(), out_of_plane_edges.end(),
                [&] (auto e1, auto e2)
                {
                    return get_max_endpoint_degree(e1) < get_max_endpoint_degree(e2);
                });
    for (sptr<net::phys_edge_t> e : out_of_plane_edges) {
        test_and_move_edge_to_bulk(e);
    }
}

inline bool
PhysicalNetwork::is_planar() {
    return lemon::checkPlanarity(planar_repr);
}

inline size_t
PhysicalNetwork::determine_edge_tsv_layer(sptr<net::phys_edge_t> e) {
    sptr<net::phys_vertex_t> src = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->src),
                            dst = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->dst);
    size_t k = 0;
    while (occupied_tsvs[src].count(k) || occupied_tsvs[dst].count(k)) k++;
    return k;
}

inline size_t
PhysicalNetwork::get_max_endpoint_degree(sptr<net::phys_edge_t> e) {
    sptr<net::phys_vertex_t> src = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->src),
                            dst = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->dst);
    const size_t sdg = get_bulk_degree(src), 
                    ddg = get_bulk_degree(dst);
    return sdg > ddg ? sdg : ddg;
}

inline size_t
PhysicalNetwork::get_bulk_degree(sptr<net::phys_vertex_t> v) {
    update_state();
    return bulk_degree_map[v];
}

}   // protean
}   // qontra
