/*
 *  author: Suhas Vittal
 *  date:   13 November 2023
 * */

#include "protean/network.h"

namespace qontra {

using namespace graph;

namespace protean {

//
//  SupportGraph implementation.
//

SupportGraph::SupportGraph(TannerGraph& tanner_graph) 
    :__SupportGraphParent(),
    vmap()
{
    // Every vertex in the tanner graph as a corresponding support graph vertex.
    for (auto tv : tanner_graph.get_vertices()) {
        sup_vertex_t* sv = new sup_vertex_t;
        sv->id = tv->id;
        sv->tgv = tv;
        if (tv->qubit_type == tanner::vertex_t::type::data) {
            sv->support = std::set<tanner::vertex_t*>{ tv };
            sv->pauli_type = sup_vertex_t::type::any;
        } else {
            sv->support = std::set<tanner::vertex_t*>(tanner_graph.get_neighbors(tv));
            sv->pauli_type = tv->qubit_type == tanner::vertex_t::type::xparity 
                                ? sup_vertex_t::type::x : sup_vertex_t::type::z;
        }
        vmap[tv] = sv;
        add_vertex(sv);
    }
    // Now, we want to build partial supports.
    uint64_t psid = 0;
    auto checks = tanner_graph.get_checks();
    for (uint i = 0; i < checks.size(); i++) {
        auto tc1 = checks[i];
        auto sc1 = get_vertex(tc1);

        std::set<sup_vertex_t*> any_unique(sc1->support);
        for (uint j = i+1; j < checks.size(); j++) {
            auto tc2 = checks[j];

            if (tc1->qubit_type != tc2->qubit_type) continue;
            auto sc2 = get_vertex(tc2);
            
            std::set<tanner::vertex_t*> shared_supp;
            std::set_intersection(
                    sc1->support.begin(), sc1->support.end(),
                    sc2->support.begin(), sc2->support.end(),
                    std::inserter(shared_supp, shared_supp.begin()));
            if (shared_supp.empty()) continue;

            sup_vertex_t* sv = new sup_vertex_t;
            sv->id = (1L << SUPPORT_VERTEX_GEN_BIT) | (psid++);
            sv->tgv = nullptr;
            sv->support = shared_supp;
            sv->pauli_type = sc1->pauli_type;
            add_vertex(sv);
            // Add edges as well.
            for (auto sd : shared_supp) {
                sup_edge_t* se = new sup_edge_t;
                se->src = (void*)sd;
                se->dst = (void*)sv;
                se->is_undirected = true;
                add_edge(se);

                any_unique.erase(sd);
            }
            sup_edge_t* se1 = new sup_edge_t;
            se1->src = (void*)sc1;
            se1->dst = (void*)sv;
            se1->is_undirected = true;

            sup_edge_t* se2 = new sup_edge_t;
            se2->src = (void*)sc2;
            se2->dst = (void*)sv;
            se2->is_undirected = true;

            add_edge(se1);
            add_edge(se2);
        }
        // If any data qubit is not in a shared support with another check,
        // add a direct edge between that qubit and this check.
        for (auto sd : any_unique) {
            sup_edge_t* se = new sup_edge_t;
            se->src = (void*)sd;
            se->dst = (void*)sc1;
            se->is_undirected = true;
            add_edge(se);
        }
    }
}

//
// Flag Proxy Network Implementation.
//

FlagProxyNetwork::FlagProxyNetwork(TannerGraph& tanner_graph)
    :__FlagProxyNetworkParent(),
    support_graph(tanner_graph),
    tanner_graph(tanner_graph),
    vmap(),
    genid(0)
{
    // Essentially, we are copying the tanner graph to here. This will be the baseline
    // connectivity.
    for (auto tv : tanner_graph.get_vertices()) {
        auto sv = support_graph.get_vertex(tv);
        auto fpv = new fp_vertex_t;
        fpv->id = tv->id;
        fpv->tanner_vertex = tv;
        fpv->support_vertex = sv;
        fpv->qubit_type = tv->qubit_type == tanner::vertex_t::type::data 
                            ? fp_vertex_t::type::data : fp_vertex_t::type::parity;
        add_vertex(fpv); vmap[sv] = fpv;
    }

    for (auto te : tanner_graph.get_edges()) {
        auto tv = (tanner::vertex_t*)te->src;
        auto tw = (tanner::vertex_t*)te->dst;
        auto fpv = get_vertex(tv);
        auto fpw = get_vertex(tw);

        if (fpv->qubit_type != fp_vertex_t::type::data) std::swap(fpv, fpw);
        // We will be using directed edges to represent the flow of data throughout
        // the network.
        auto fpe = new fp_edge_t;
        fpe->src = (void*)fpv;
        fpe->dst = (void*)fpw;
        fpe->is_undirected = false;
        add_edge(fpe);
    }
}

void
FlagProxyNetwork::optimize() {
    // First, we need to require that the graph meets connectivity constraints.
}

void
FlagProxyNetwork::reduce_connectivity_of_data(fp_vertex_t* fpv) {
    auto tv = fpv->tanner_vertex;
    // Find other data qubits that are violating the connectivity requirements
    // that are connected to one of the parity checks of fpv.
    tanner::vertex_t* check_w_most_incidence = nullptr;
    uint max_incidence_cnt = 0;
    std::set<fp_vertex_t*> max_incidence_violators;
    for (auto tc : tanner_graph.get_neighbors(tv)) {
        if (is_flagged_map[tc].count(tv)) continue;

        uint cnt = 0;
        std::set<fp_vertex_t*> violators;
        for (auto tw : tanner_graph.get_neighbors(tc)) {
            if (tv == tw) continue;
            fp_vertex_t* fpw = get_vertex(tw);
            // We will mark fpw as a violator if it violates connectivity
            // requirements AND it has a direct connection to the parity qubit.
            if (get_degree(fpw) > config.min_connectivity
                && !is_flagged_map[tc].count(tw))
            {
                violators.insert(fpw);
                cnt++;
            }
        }
        if (cnt > max_incidence_cnt) {
            check_w_most_incidence = tc;
            max_incidence_cnt = cnt;
            max_incidence_violators = violators;
        }
    }
    // If there are no violators that share a check with the data qubit,
    // it is by itself. Just introduce proxies.
    if (max_incidence_violators.empty()) {
        reduce_connectivity_by_adding_proxy(fpv);
        return;
    }
    // Otherwise, try and make a flag.
    auto flag = new fp_vertex_t;
    flag->id = (1L << FLAG_PROXY_VERTEX_GEN_BIT) | (genid++);
    flag->qubit_type = fp_vertex_t::type::flag;
    add_vertex(flag);

    violators.insert(fpv);
    for (auto fpx : violators) {
        auto e = new fp_edge_t;
        e->src = (void*)fpx;
        e->dst = (void*)flag;
        e->is_undirected = true;
        add_edge(e);
    }
}

void
FlagProxyNetwork::reduce_connectivity_of_parity(fp_vertex_t* fpv) {
    auto sv = fpv->support_vertex;
}

void
FlagProxyNetwork::reduce_connectivity_by_adding_proxies(fp_vertex_t* fpv) {
    add_proxies_to(fpv, get_indegree(fpv) > get_outdegree(fpv));
    add_proxies_to(fpv, get_indegree(fpv) <= get_outdegree(fpv));
}

void
FlagProxyNetwork::add_proxies_to(fp_vertex_t* fpv, bool on_incoming) {
    int neg_slack = get_inoutdegree(fpv) - config.min_connectivity;
    if (neg_slack <= 0) return; // Nothing to be done.
    // Create proxy:
    auto curr_proxy = make_proxy(fpv, on_incoming);
    neg_slack++;    // The slack increases as we added a proxy edge.
    // Add edges to proxy.
    auto neighbors = on_incoming ? get_incoming(fpv) : get_outgoing(fpv);
    for (auto fpw : neighbors) {
        auto fpe_old = on_incoming ? get_edge(fpw, fpv) : get_edge(fpv, fpw);
        // Delete this edge.
        delete_edge(fpe_old);
        // Create new edge with the proxy.
        auto fpe_new = new fp_edge_t;
        fpe_new->src = (void*)(on_incoming ? fpw : curr_proxy);
        fpe_new->dst = (void*)(on_incoming ? curr_proxy : fpw);
        fpe_new->is_undirected = false;
        add_edge(fpe_new);
        neg_slack--;
        if (neg_slack == 0) break;  // We are done.
        if (get_inoutdegree(curr_proxy) == config.min_connectivity) {
            // Get a new proxy.
            curr_proxy = make_proxy(fpv, on_incoming);
            neg_slack++;
        }
    }
    // If the current proxy has done nothing, delete it.
    if (get_inoutdegree(curr_proxy) == 1) {
        delete_vertex(curr_proxy);
    }
}

fp_vertex_t*
FlagProxyNetwork::make_proxy(fp_vertex_t* fpv, bool incoming) {
    auto proxy = new fp_vertex_t;
    proxy->id = genid++;
    proxy->qubit_type = fp_vertex_t::type::proxy;
    add_vertex(proxy);
    // Add edge between fpv and proxy.
    auto e = new fp_edge_t;
    e->src = (void*)(incoming ? proxy : fpv);
    e->dst = (void*)(incoming ? fpv : proxy);
    return proxy;
}

}   // protean
}   // qontra
