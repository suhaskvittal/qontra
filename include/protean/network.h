/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#ifndef PROTEAN_NETWORK_h
#define PROTEAN_NETWORK_h

#include <defs/bijective_map.h>
#include <defs/two_level_map.h>
#include <graph/tanner_graph.h>

#include <lemon/list_graph.h>

namespace qontra {
namespace protean {

// We will have two types of networks:
//  (1) A "raw" network, which has the following rules:
//      (i) a qubit can only be used for "one" role (i.e. parity qubit for a single check, a flag, etc.)
//  (2) A physical network, which is the actual coupling graph. Here, qubits can have multiple roles
//      at different points in time.

namespace net {

struct raw_vertex_t : graph::base::vertex_t {
    enum class type { data, xparity, zparity, flag, proxy };
    type qubit_type;
};

struct raw_edge_t : graph::base::edge_t {};

struct phys_vertex_t : graph::base::vertex_t {
    std::set<sptr<raw_vertex_t>> role_set;
    BijectiveMap<size_t, sptr<raw_vertex_t>> cycle_role_map;

    void add_role(sptr<raw_vertex_t>, size_t cycle);
    void push_back_role(sptr<raw_vertex_t>);
    bool has_role_of_type(raw_vertex_t::type);
private:
    std::set<raw_vertex_t::type> role_type_set;
};

struct phys_edge_t : graph::base::edge_t {
    bool is_out_of_plane;
    size_t tsv_layer;
};

}   // net

class RawNetwork : public graph::Graph<net::raw_vertex_t, net::raw_edge_t> {
public:
    RawNetwork(graph::TannerGraph&);

    // Input: data qubit, data qubit, parity qubit
    // Returns: flag qubit.
    sptr<net::raw_vertex_t> add_flag(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // Input: xyz qubit, xyz qubit (xyz meaning can be any type).
    // Returns: proxy qubit. The proxy qubit splits the edge between the two input qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // This does the same as the above command, but just takes an edge instead of two qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_edge_t>);

    // Input: xyz qubit, proxy qubit, a reference to a vector which will be populated with the walk.
    // Returns: non-proxy qubit obtained by walking through the proxy_indirection_map.
    sptr<net::raw_vertex_t> proxy_walk(sptr<net::raw_vertex_t>,
                                            sptr<net::raw_vertex_t>,
                                            std::vector<sptr<net::raw_vertex_t>>&);
    
    // Tanner graph tracking structures:
    BijectiveMap<sptr<graph::tanner::vertex_t>, sptr<net::raw_vertex_t>> v_tanner_raw_map;

    // Flag tracking structures:
    //
    // parity qubit --> list of flag qubits used in syndrome extraction of check.
    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>> 
        flag_ownership_map;
    // parity qubit --> data qubit --> flag qubit used in syndrome extraction of check.
    TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        flag_assignment_map;
    // proxy qubit --> xyz qubit --> qubit xyz was linked to before
    TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        proxy_indirection_map;

    TannerGraph tanner_graph;
private:
    uint64_t id_ctr = 0;
};

class PhysicalNetwork : public graph::Graph<net::phys_vertex_t, net::phys_edge_t> {
public:
    PhysicalNetwork(graph::TannerGraph&, uint max_connectivity);
    // All of the graph modification functions need to be updated to handle planarity.
    // Note that deletes only free up spots on the processor bulk. Only adding edges
    // requires checking for planarity.
    bool add_vertex(sptr<net::phys_vertex_t>) override;
    bool add_edge(sptr<net::phys_edge_t>) override;
    void delete_vertex(sptr<net::phys_vertex_t>) override;
    void delete_edge(sptr<net::phys_edge_t>) override;

    // Returns true if the graph remains planar.
    bool test_and_move_edge_to_bulk(sptr<net::phys_edge_t>, bool is_new_edge=false); 

    void make_flags(bool for_x_parity);
protected:
    bool update_state(void) override;
private:
    void reallocate_edges(void); // Try and place out-of-plane edges into the processor bulk.
    bool is_planar(void);
    size_t determine_edge_tsv_layer(sptr<net::phys_edge_t>);

    size_t get_max_endpoint_degree(sptr<net::phys_edge_t>);
    size_t get_bulk_degree(sptr<net::phys_vertex_t>);

    // raw_connection_network contains all roles in the network, from proxy to flag to data (etc.).
    // Each phys_vertex_t corresponds to at least one raw_vertex_t (if not more).
    RawNetwork raw_connection_network;
    std::map<net::raw_vertex_t, net::phys_vertex_t> role_to_phys;
    // planar_representation will contain all edges on the bulk of the processor (not
    // edges going through a TSV). If an edge can be added to the planar_representation
    // without violating planarity, it goes on the bulk.
    lemon::ListGraph planar_repr;
    BijectiveMap<sptr<net::phys_vertex_t>, lemon::ListGraph::Node>  v_phys_lemon_map;
    BijectiveMap<sptr<net::phys_edge_t>, lemon::ListGraph::Edge>    e_phys_lemon_map;
    // Other tracking structures:
    //
    // Tracks the heights of TSV edges for each vertex.
    std::map<sptr<net::phys_vertex_t>, std::set<size_t>> occupied_tsvs;
    // Tracks the degree of each vertex in the processor bulk.
    std::map<sptr<net::phys_vertex_t>, size_t> bulk_degree_map;

    const uint max_connectivity;
};

}   // protean
}   // qontra

#include "network.inl"

#endif  // PROTEAN_NETWORK_h
