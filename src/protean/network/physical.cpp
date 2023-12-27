/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include "protean/network.h"

namespace qontra {
namespace protean {

using namespace graph;
using namespace net;

PhysicalNetwork(TannerGraph& tgr, uint max_conn)
    :Graph(),
    raw_connection_network(tgr),
    role_to_phys(),
    // Planarity tracking:
    planar_repr(),
    v_phys_to_lemon(),
    v_lemon_to_phys(),
    e_phys_to_lemon(),
    e_lemon_to_phys(),
    // Other tracking:
    occupied_tsvs(),
    bulk_degree_map(),
    // Other variables:
    max_connectivity(max_conn)
{
    // Create a corresponding physical qubit for every vertex in raw_connection_network.
    for (sptr<raw_vertex_t> rv : raw_connection_network.get_vertices()) {
        sptr<phys_vertex_t> pv = std::make_shared<>();
        pv->id = rv->id;
        pv->push_back_role(rv);
        role_to_phys[rv] = pv;
        add_vertex(pv);
    }
    // And the corresponding edges as well.
    for (sptr<raw_edge_t> re : raw_connection_network.get_edges()) {
        sptr<raw_vertex_t> rsrc = std::reinterpret_pointer_cast<raw_vertex_t>(re->src),
                            rdst = std::reinterpret_pointer_cast<raw_vertex_t>(re->dst);
        sptr<phys_edge_t> pe = std::make_shared<>();
        pe->src = role_to_phys[rsrc];
        pe->dst = role_to_phys[rdst];
        pe->is_undirected = true;
        add_edge(pe);
    }
}

void
PhysicalNetwork::make_flags() {

}

}   // protean
}   // qontra
