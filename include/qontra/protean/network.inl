/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include <qontra/graph/algorithms/planarity.h>

#include <algorithm>

namespace qontra {

namespace graph {

template <> inline std::string
print_v<protean::net::raw_vertex_t>(sptr<protean::net::raw_vertex_t> v) {
    using namespace protean;
    using namespace net;
    if (v == nullptr) return "<nil>";

    std::string s;
    if (v->qubit_type == raw_vertex_t::type::data) {
        s += "d";
    } else if (v->qubit_type == raw_vertex_t::type::xparity) {
        s += "x";
    } else if (v->qubit_type == raw_vertex_t::type::zparity) {
        s += "z";
    } else if (v->qubit_type == raw_vertex_t::type::flag) {
        s += "f";
    } else {
        s += "pr";
    }
    s += std::to_string(v->id);
    return s;
}

} // graph

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
    role_set.insert(r);
    role_type_map[r->qubit_type]++;
}

inline void
phys_vertex_t::delete_role(sptr<raw_vertex_t> r) {
    role_set.erase(r);
    cycle_role_map.erase(r);
    if ((--role_type_map[r->qubit_type]) == 0) role_type_map.erase(r->qubit_type);
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

inline void
phys_vertex_t::clear_roles() {
    cycle_role_map.clear();
    role_set.clear();
}

inline bool
phys_vertex_t::has_role_of_type(raw_vertex_t::type t) {
    return role_type_map.count(t);
}

inline bool
phys_edge_t::is_out_of_plane() {
    return tsv_layer > 0;
}

}   // net

inline sptr<net::raw_vertex_t>
RawNetwork::make_vertex() {
    sptr<net::raw_vertex_t> v = std::make_shared<net::raw_vertex_t>();
    v->id = id_ctr++;
    return v;
}

inline sptr<net::raw_vertex_t>
RawNetwork::add_proxy(sptr<net::raw_vertex_t> q1, sptr<net::raw_vertex_t> q2) {
    sptr<net::raw_edge_t> e = get_edge(q1, q2);
    return add_proxy(e);
}

inline sptr<net::raw_vertex_t>
RawNetwork::are_in_same_support(sptr<net::raw_vertex_t> rx, sptr<net::raw_vertex_t> ry) {
    if (!enable_memoization ||
            (!same_support_memo_map.count(rx) || !same_support_memo_map[rx].count(ry)))
    {
        auto rpq = are_in_same_support({rx, ry});
        same_support_memo_map[rx][ry] = rpq;
        same_support_memo_map[ry][rx] = rpq;
    }
    return same_support_memo_map[rx][ry];
}

inline void
RawNetwork::reset_memoization() {
    proxy_memo_map.clear();
    support_memo_map.clear();
    same_support_memo_map.clear();
}

inline void
RawNetwork::disable_memoization() {
    reset_memoization();
    enable_memoization = false;
}

inline bool
ProcessorLayer::add_edge(sptr<net::phys_edge_t> e) {
    if (!Graph::add_edge(e)) return false;
    // Otherwise, the edge has been added. Check planarity.
    if (!is_planar()) {
        // Delete the edge.
        delete_edge(e);
        return false;
    }
    return true;
}

inline size_t
ProcessorLayer::get_max_endpoint_degree(sptr<net::phys_edge_t> e) {
    sptr<net::phys_vertex_t> v = get_source(e),
                             w = get_target(e);
    return std::max(get_degree(v), get_degree(w));
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

inline sptr<net::phys_vertex_t>
PhysicalNetwork::get_physical_qubit_for(sptr<net::raw_vertex_t> r) {
    return role_to_phys[r];
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
    for (size_t k = 0; k < get_thickness(); k++) {
        auto& layer = processor_layers[k];
        if (layer.is_planar() && layer.add_edge(e)) {
            e->tsv_layer = k;
            return true;
        }
    }
    // If none can be found, then make a new layer.
    std::cout << "[ debug ] creating new processor layer\n";
    ProcessorLayer& new_layer = push_back_new_processor_layer();
    if (new_layer.add_edge(e)) {
        e->tsv_layer = get_thickness()-1;
        return true;
    } else {
        sptr<net::phys_vertex_t> src = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->src),
                                 dst = std::reinterpret_pointer_cast<net::phys_vertex_t>(e->dst);
        std::cout << "[ debug ] failed to add edge " << graph::print_e<net::phys_vertex_t>(e)
                        << ", reasons:\n"
                        << "\tlayer does not contain endpoints: "
                        << graph::print_v(src) << "(" << new_layer.contains(src)
                        << "), " << graph::print_v(dst) << "(" << new_layer.contains(dst) << ")\n";
        Graph::delete_edge(e);
        return false;
    }
}

inline void
PhysicalNetwork::delete_vertex(sptr<net::phys_vertex_t> v) {
    Graph::delete_vertex(v);
    for (auto& layer : processor_layers) {
        layer.delete_vertex(v);
    } 
}

inline void
PhysicalNetwork::delete_edge(sptr<net::phys_edge_t> e) {
    Graph::delete_edge(e);
    processor_layers[e->tsv_layer].delete_edge(e);
}

inline sptr<net::phys_vertex_t>
PhysicalNetwork::make_vertex() {
    sptr<net::phys_vertex_t> v = std::make_shared<net::phys_vertex_t>();
    v->id = id_ctr++;
    return v;
}

inline bool
PhysicalNetwork::test_and_move_edge_down(sptr<net::phys_edge_t> e) {
    if (!e->is_out_of_plane()) return false;
    auto& curr_layer = processor_layers[e->tsv_layer];
    auto& next_layer = processor_layers[e->tsv_layer-1];

    if (next_layer.add_edge(e)) {
        // Move was successful.
        curr_layer.delete_edge(e);
        e->tsv_layer--;
        return true;
    } else {
        return false;
    }
}

inline size_t
PhysicalNetwork::get_thickness() {
    return processor_layers.size();
}

inline size_t
PhysicalNetwork::get_bulk_degree(sptr<net::phys_vertex_t> v) {
    return processor_layers[0].get_degree(v);
}

inline bool
PhysicalNetwork::relabel_qubits() {
    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<net::phys_vertex_t> pv = vertices[i];
        if (pv->id == i) continue;
        // Otherwise, relabel the vertex.
        manual_update_id(pv, pv->id, i);
        for (auto& layer : processor_layers) {
            layer.manual_update_id(pv, pv->id, i);
        }
        pv->id = i;
    }
    id_ctr = vertices.size();
    return true;
}

inline ProcessorLayer&
PhysicalNetwork::push_back_new_processor_layer() {
    ProcessorLayer new_layer;
    for (sptr<net::phys_vertex_t> v : get_vertices()) {
        new_layer.add_vertex(v);
    }
    processor_layers.push_back(new_layer);
    return processor_layers.back();
}

inline void
PhysicalNetwork::consume(sptr<net::phys_vertex_t> v, sptr<net::phys_vertex_t> w) {
    for (auto r : w->role_set) role_to_phys[r] = v;
    v->consume(w);
}

inline bool
PhysicalNetwork::are_in_same_support(sptr<net::phys_vertex_t> v, sptr<net::phys_vertex_t> w) {
    for (auto rv : v->role_set) {
        for (auto rw : w->role_set) {
            if (raw_connection_network.are_in_same_support(rv, rw) != nullptr) return true;
        }
    }
    return false;
}

}   // protean
}   // qontra
