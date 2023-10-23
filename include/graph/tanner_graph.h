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

const int VERTEX_ID_AUTOGEN_BIT = 63;

struct vertex_t : graph::base::vertex_t {
    enum class Type { data, xparity, zparity };
    Type qubit_type;
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
        x_parity_checks(),
        z_parity_checks(),
        x_obs_list(),
        z_obs_list()
    {}

    TannerGraph(const TannerGraph& other)
        :__TannerGraphParent(other),
        data_qubits(other.data_qubits),
        x_parity_checks(other.x_parity_checks),
        z_parity_checks(other.z_parity_checks),
        x_obs_list(other.x_obs_list),
        z_obs_list(other.z_obs_list)
    {}

    bool add_vertex(tanner::vertex_t* v) override {
        if (!__TannerGraphParent::add_vertex(v))  return false;
        if (v->qubit_type == tanner::vertex_t::Type::data)      data_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::Type::xparity)   x_parity_checks.push_back(v);
        if (v->qubit_type == tanner::vertex_t::Type::zparity)   z_parity_checks.push_back(v);
        return true;
    }

    bool add_edge(tanner::edge_t* e) override {
        auto v = (tanner::vertex_t*)e->src;
        auto w = (tanner::vertex_t*)e->dst;
        // Make sure the edge preserves the bipartite property.
        bool src_is_parity = v->qubit_type != tanner::vertex_t::Type::data;
        bool dst_is_parity = v->qubit_type != tanner::vertex_t::Type::data;
        if (src_is_parity && dst_is_parity) return false;
        return __TannerGraphParent::add_edge(e);
    }

    void delete_vertex(tanner::vertex_t* v) override {
        std::vector<tanner::vertex_t*>* cat;
        if (v->qubit_type == tanner::vertex_t::Type::data)      cat = &data_qubits;
        if (v->qubit_type == tanner::vertex_t::Type::xparity)   cat = &x_parity_qubits;
        if (v->qubit_type == tanner::vertex_t::Type::zparity)   cat = &z_parity_qubits;
        for (auto it = cat->begin(); it != cat->end();) {
            if (*it == v)   it = cat->erase(it);
            else            it++;
        }
        __TannerGraphParent::delete_vertex(v);
    }

    typedef std::vector<tanner::vertex_t*>  obs_t;

    std::vector<obs_t>  x_obs_list;
    std::vector<obs_t>  z_obs_list;
private:
    std::vector<tanner::vertex_t*>  data_qubits;
    std::vector<tanner::vertex_t*>  x_parity_checks;
    std::vector<tanner::vertex_t*>  z_parity_checks;
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
