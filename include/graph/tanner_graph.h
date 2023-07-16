/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#ifndef TANNER_GRAPH_h
#define TANNER_GRAPH_h

#include "defs.h"
#include "graph/graph.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace qontra {
namespace graph {

namespace tanner {

struct vertex_t : graph::base::vertex_t {
    uint8_t qubit_type; // 00 = Data, 01 = Gauge, 11 = X parity, 10 = Z parity

    const static uint8_t DATA =     0x0;
    const static uint8_t GAUGE =    0x1;
    const static uint8_t XPARITY =  0x2;
    const static uint8_t ZPARITY =  0x3;
};

struct edge_t : graph::base::edge_t {
};

}   // tanner

#define __TannerGraphParent graph::Graph<tanner::vertex_t, tanner::edge_t>

class TannerGraph : public __TannerGraphParent {
public:
    TannerGraph(void)
        :__TannerGraphParent(), 
        data_qubits(),
        gauge_qubits(),
        x_parity_checks(),
        z_parity_checks(), 
        induced_gauge_index(0)
    {}

    TannerGraph(const TannerGraph& other)
        :__TannerGraphParent(other),
        data_qubits(other.data_qubits),
        gauge_qubits(other.gauge_qubits),
        x_parity_checks(other.x_parity_checks),
        z_parity_checks(other.z_parity_checks),
        induced_gauge_index(other.induced_gauge_index)
    {}

    bool add_vertex(tanner::vertex_t* v) override {
        if (!Graph::add_vertex(v))  return false;
        if (v->qubit_type == tanner::vertex_t::DATA)     data_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::GAUGE)    gauge_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::XPARITY)  x_parity_checks.push_back(v);
        if (v->qubit_type == tanner::vertex_t::ZPARITY)  z_parity_checks.push_back(v);
        return true;
    }

    bool add_edge(tanner::edge_t* e) override {
        auto v = (tanner::vertex_t*)e->src;
        auto w = (tanner::vertex_t*)e->dst;
        // Make sure the edge preserves the bipartite property.
        if ((v->qubit_type > 0) == (w->qubit_type > 0)) return false;
        bool res = Graph::add_edge(e);
        if (res) {
            // Sort adjacency lists.
            std::sort(adjacency_lists[v].begin(), adjacency_lists[v].end());
            if (e->is_undirected)   std::sort(adjacency_lists[w].begin(), adjacency_lists[w].end());
        }
        return res;
    }

    void delete_vertex(tanner::vertex_t* v) override {
        std::vector<tanner::vertex_t*>* cat;
        if (v->qubit_type == tanner::vertex_t::DATA)     cat = &data_qubits;
        if (v->qubit_type == tanner::vertex_t::GAUGE)    cat = &gauge_qubits;
        if (v->qubit_type == tanner::vertex_t::XPARITY)  cat = &x_parity_checks;
        if (v->qubit_type == tanner::vertex_t::ZPARITY)  cat = &z_parity_checks;
        for (auto it = cat->begin(); it != cat->end();) {
            if (*it == v)   it = cat->erase(it);
            else            it++;
        }
        __TannerGraphParent::delete_vertex(v);
    }

    std::vector<tanner::vertex_t*>  get_predecessors(tanner::vertex_t*); 
                                                                // Gets the set of predecessor 
                                                                // for a vertex. A vertex X 
                                                                // precedes a vertex Y if
                                                                // adj(X) is a subset of adj(Y).
    tanner::vertex_t*       induce_predecessor(tanner::vertex_t*, tanner::vertex_t*);  
                                                                // An induced predecessor of 
                                                                // two vertices X and Y is a 
                                                                // vertex Z such that adj(Z) 
                                                                // is the intersection of 
                                                                // adj(X) and adj(Y). The induced
                                                                // predecessor is added to the
                                                                // graph. This function fails
                                                                // if either arg precedes the
                                                                // other.
    bool                    has_copy_in_gauges(const std::vector<tanner::vertex_t*>& adj);
    
    std::vector<tanner::vertex_t*>  get_vertices_by_type(uint8_t t) {
        const std::vector<tanner::vertex_t*> arr[] = 
                { data_qubits, gauge_qubits, x_parity_checks, z_parity_checks };
        return arr[t];
    }
private:
    std::vector<tanner::vertex_t*>  data_qubits;
    std::vector<tanner::vertex_t*>  gauge_qubits;
    std::vector<tanner::vertex_t*>  x_parity_checks;
    std::vector<tanner::vertex_t*>  z_parity_checks;

    uint induced_gauge_index;
    const static uint INDUCED_GAUGE_INDEX_FLAG = 1 << 24;
};

namespace io {

// Function for io callback.
//
// In a tanner graph description file, each line should be of the form:
//          <X/Z><check-id>,<data-qubit-1>,<data-qubit-2>,...

void    update_tanner_graph(TannerGraph&, std::string); // Callback for io function.

}   // io

}   // graph
}   // qontra

#endif  // TANNER_GRAPH_h
