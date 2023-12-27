/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include <lemon/planarity.h>

namespace qontra {
namespace protean {

namespace net {

inline void
phys_vertex_t::push_back_role(sptr<raw_vertex_t> r) {
    const uint cycle = role_set.size();
    role_set.insert(r);
    cycle_to_role[cycle] = r;
    role_to_cycle[r] = cycle;
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
        phys_to_lemon[v] = lemon_node;
        v_lemon_to_phys[lemon_node] = v;
        return true;
    } else {
        return false;
    }
}

inline bool
PhysicalNetwork::add_edge(sptr<net::phys_edge_t> e) {
    if (graph::Graph::add_edge(e)) {
        sptr<net::phys_vertex_t> src = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->src),
                                dst = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->dst);
        auto lemon_src = v_phys_to_lemon[src],
                lemon_dst = v_phys_to_lemon[dst];
        auto lemon_edge = planar_repr.addEdge(lemon_src, lemon_dst);
        if (!is_planar()) {
            planar_repr.erase(lemon_edge);
            e->is_out_of_plane = true;

            const size_t k = determine_edge_tsv_layer(e);
            e->tsv_layer = k;
            occupied_tsvs[src].insert(k);
            occupied_tsvs[dst].insert(k);
        } else {
            e_phys_to_lemon[e] = lemon_edge;
            e_lemon_to_phys[lemon_edge] = e;
            e->is_out_of_plane = false;
        }
        return true;
    } else {
        return false;
    }
}

inline void
PhysicalNetwork::delete_vertex(sptr<net::phys_vertex_t> v) {
    auto lemon_node = v_phys_to_lemon[v];
    // First, erase all traces of v in the tracking structures.
    v_phys_to_lemon.erase(v);
    v_lemon_to_phys.erase(lemon_node);
    for (sptr<net::phys_edge_t> e : get_neighbors(v)) {
        if (e->is_out_of_plane) continue;
        auto lemon_edge = e_phys_to_lemon[e];

        e_phys_to_lemon.erase(e);
        e_lemon_to_phys.erase(lemon_edge);

        planar_repr.erase(lemon_edge);
    }
    // Then delete lemon_node and v.
    planar_repr.erase(lemon_node);
    graph::Graph::delete_vertex(v);
}

inline void
PhysicalNetwork::delete_edge(sptr<net::phys_edge_t> e) {
    if (!e->is_out_of_plane) {
        auto lemon_edge = e_phys_to_lemon[e];
        e_phys_to_lemon.erase(e);
        e_lemon_to_phys.erase(lemon_edge);
        planar_repr.erase(lemon_edge);
    }
    graph::Graph::delete_edge(e);
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

}   // protean
}   // qontra
