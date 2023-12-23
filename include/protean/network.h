/*
 *  author: Suhas Vittal
 *  date:   10 November 2023
 * */

#ifndef PROTEAN_NETWORK_h
#define PROTEAN_NETWORK_h

#include "graph/graph.h"
#include "graph/tanner_graph.h"

#include <algorithm>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

//
// Support Graph
//
// The SupportGraph encodes relationships between data qubits and "support pools",
// or sets of qubits. This is an extension of a TannerGraph, where two stabilizers
// (of the same type) sharing a suppport will have a child X which contains the 
// intersection of their supports. Any data qubits are connected to this intersection.
//

const int SUPPORT_VERTEX_GEN_BIT = 60;

struct sup_vertex_t : graph::base::vertex_t {
    enum class type { any, x, z };

    tanner::vertex_t*                   tgv;
    std::set<graph::tanner::vertex_t*>  support;
    type pauli_type;
};

struct sup_edge_t : graph::base::edge_t {
};

#define __SupportGraphParent graph::Graph<sup_vertex_t, sup_edge_t>
class SupportGraph : public __SupportGraphParent {
public:
    SupportGraph(TannerGraph&);
    SupportGraph(const SupportGraph& other)
        :__SupportGraphParent(other),
        vmap(other.vmap)
    {}

    sup_vertex_t* get_vertex(graph::tanner::vertex_t* tv) {
        if (vmap.count(tv)) return vmap[tv];
        else                return nullptr;
    }
private:
    std::map<graph::tanner::vertex_t*, sup_vertex_t*> vmap;
};

//
// Flag Proxy Network.
//

const int FLAG_PROXY_VERTEX_GEN_BIT = 48;

struct fp_vertex_t : graph::base::vertex_t {
    enum class type { data, flag, proxy, parity };

    graph::tanner::vertex_t*    tanner_vertex;
    sup_vertex_t*               support_vertex;

    type qubit_type;
};

struct fp_edge_t : graph::base::edge_t {
};

#define __FlagProxyNetworkParent    graph::Graph<fp_vertex_t, fp_edge_t>
class FlagProxyNetwork : public __FlagProxyNetworkParent {
public:
    FlagProxyNetwork(TannerGraph&);
    FlagProxyNetwork(const FlagProxyNetwork& other)
        :__FlagProxyNetworKParent(other),
        support_graph(other.support_graph),
        tanner_graph(other.tanner_graph),
        vmap(other.vmap),
        proxyid(other.proxyid)
    {}

    struct {
        // Constraints
        int min_connectivity = 4;
        // Optimization setting.
        enum class opt { depth, conn, qubits };
        opt target = opt::conn;
    } config;

    void optimize(void);

    fp_vertex_t* get_vertex(graph::tanner::vertex_t* tv) {
        return get_vertex(support_graph.get_vertex(tv));
    }

    fp_vertex_t* get_vertex(sup_vertex_t* sv) {
        if (vmap.count(sv)) return vmap[sv];
        return nullptr;
    }
private:
    void add_partial_support_to(

    void reduce_connectivity(fp_vertex_t*);

    void reduce_connectivity_of_data(fp_vertex_t*);
    void reduce_connectivity_of_parity(fp_vertex_t*);
    void reduce_connectivity_of_flag(fp_vertex_t*);

    void reduce_connectivity_by_adding_proxies(fp_vertex_t*);

    void add_proxies_to(fp_vertex_t*, bool on_incoming);

    fp_vertex_t* make_proxy(fp_vertex_t*, bool incoming);

    SupportGraph        support_graph;
    graph::TannerGraph  tanner_graph;

    std::map<sup_vertex_t*, fp_vertex_t*> vmap;

    std::map<tanner::vertex_t*, std::set<tanner::vertex_t*>> is_flagged_map;

    uint64_t genid;
};

}   // protean
}   // qontra

#endif  // PROTEAN_NETWORK_h
