/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include <defs/filesystem.h>
#include <graph/algorithms/planarity.h>

#include <algorithm>


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
ProcessorLayer::add_edge(sptr<net::phys_edge_t> e) {
    if (!Graph::add_edge(e)) return false;
    // Otherwise, the edge has been added. Check planarity.
    if (!is_planar()) {
        // Delete the edge.
        delete_edge(e);
        return e;
    }
    return true;
}

inline bool
ProcessorLayer::is_planar() {
    update_state();
    return _is_planar;
}

inline bool
ProcessorLayer::update_state() {
    if (!Graph::update_state()) return false;
    // Update planarity.
    _is_planar = lr_planarity(this);
    return true;
}

inline RawNetwork
PhysicalNetwork::get_raw_connection_network() {
    return raw_connection_network;
}

inline bool
PhysicalNetwork::add_vertex(sptr<net::phys_vertex_t> v) {
    if (!Graph::add_vertex(v)) return false;
    // Add v to each processor layer.
    for (ProcessorLayer& pl : processor_layers) pl.add_vertex(v);
    return true;
}

inline bool
PhysicalNetwork::add_edge(sptr<net::phys_edge_t> e) {
    if (!Graph::add_edge(e)) return false;
    // Here, we will attempt to find a processor layer that will house e.
    for (size_t k = 0; k < processor_layers.size(); k++) {
        if (processor_layers[k].add_edge(e)) return true;
    }
    // If none can be found, then make a new layer.
    ProcessorLayer& new_layer = push_back_new_processor_layer();
    new_layer.add_edge(e);
    return true;
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
    return processor_layers[0].get_degree(v);
}

inline void
write_network_to_folder(std::string output_folder, PhysicalNetwork& network) {
    // Make folder if it does not exist.
    safe_create_directory(output_folder);

    std::string schedule_file = output_folder + "/round.asm";
    std::string coupling_file = output_folder + "/coupling_graph.txt";
    std::string role_file = output_folder + "/roles.txt";
    std::string tanner_graph_file = output_folder + "/tanner_graph.txt";
    std::string flag_assign_file = output_folder + "/flag_assignment.txt";

    write_schedule_file(schedule_file, network);
    write_coupling_file(coupling_file, network);
    write_role_file(role_file, network);
    write_tanner_graph_file(tanner_graph_file, network);
    write_flag_assignment_file(flag_assign_file, network);
}

}   // protean
}   // qontra
