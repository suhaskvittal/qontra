/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#ifndef PROTEAN_TANNER_GRAPH_h
#define PROTEAN_TANNER_GRAPH_h

#include "defs.h"
#include "graph/graph.h"

#include <algorithm>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

namespace tanner {

struct vertex_t : base::vertex_t {
    uint8_t qubit_type; // 00 = Data, 01 = Gauge, 11 = X parity, 10 = Z parity

    const uint8_t TANNER_DATA =     0x0;
    const uint8_t TANNER_GAUGE =    0x1;
    const uint8_t TANNER_XPARITY =  0x2;
    const uint8_t TANNER_ZPARITY =  0x3;
};

struct edge_t : base::edge_t<vertex_t> {
};

}   // tanner

class TannerGraph : Graph<tanner::vertex_t, tanner::edge_t> {
public:
    TannerGraph(void)
        :Graph(), data_qubits(), gauge_qubits(), x_parity_checks(), z_parity_checks(), 
        induced_gauge_index(0)
    {}

    TannerGraph(const TannerGraph& other)
        :Graph(other),
        data_qubits(other.data_qubits),
        gauage_qubits(other.gauge_qubits),
        x_parity_checks(other.x_parity_checks),
        z_parity_checks(other.z_parity_checks),
        induced_gauge_index(other.induced_gauge_index)
    {}

    bool add_vertex(tanner::vertex_t* v) override {
        if (v->qubit_type == tanner::vertex_t::TANNER_DATA)     data_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::TANNER_GAUGE)    gauge_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::TANNER_XPARITY)  x_parity_checks.push_back(v);
        if (v->qubit_type == tanner::vertex_t::TANNER_ZPARITY)  z_parity_checks.push_back(v);
        return Graph::add_vertex(v);
    }

    bool add_edge(tanner::edge_t* e) override {
        auto v = e->src;
        auto w = e->dst;
        // Make sure the edge preserves the bipartite property.
        if ((v->qubit_type > 0) == (w->qubit_type > 0)) return false;
        else                                            return Graph::add_edge(e);
    }

    void delete_vertex(tanner::vertex_t* v) override {
        std::vector<tanner::vertex_t*>& cat;
        if (v->qubit_type == tanner::vertex_t::TANNER_DATA)     cat = data_qubits;
        if (v->qubit_type == tanner::vertex_t::TANNER_GAUGE)    cat = gauge_qubits;
        if (v->qubit_type == tanner::vertex_t::TANNER_XPARITY)  cat = x_parity_checks;
        if (v->qubit_type == tanner::vertex_t::TANNER_ZPARITY)  cat = z_parity_checks;
        for (auto it = cat.begin(); it != cat.end();) {
            if (*it == v)   it = cat.erase(it);
            else            it++;
        }
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
                                                                // graph.

    
    std::vector<tanner::vertex_t*>  get_vertices_by_type(uint8_t t) {
        const std::vector<tanner::vertex_t*> arr[] = 
                { data_qubits, gauge_qubits, z_parity_checks, x_parity_checks };
        return arr[t];
    }
private:
    std::vector<tanner::vertex_t*>  data_qubits;
    std::vector<tanner::vertex_t*>  gauge_qubits;
    std::vector<tanner::vertex_t*>  x_parity_checks;
    std::vector<tanner::vertex_t*>  z_parity_checks;

    uint induced_gauge_index;
    const uint INDUCED_GAUGE_INDEX_FLAG = 1 << 24;
};

}   // protean
}   // qontra

#endif  // PROTEAN_TANNER_GRAPH_h
