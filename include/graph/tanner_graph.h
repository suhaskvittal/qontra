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

const uint64_t VERTEX_ID_NUMBER_MASK = (1L<<32)-1;

// The upper three bits of the tanner vertex ID are reserved
// for identification of the vertex. This is used to avoid
// conflicts when generating a tanner graph from a specification
// file.
//
// If a user is building their tanner graph from scratch, we recommend
// they use these constants when constructing ids.
const uint64_t VERTEX_ID_DATA_FLAG = (1L<<61);
const uint64_t VERTEX_ID_XPARITY_FLAG = (1L<<62);
const uint64_t VERTEX_ID_ZPARITY_FLAG = (1L<<63);

struct vertex_t : graph::base::vertex_t {
    enum class Type { data, xparity, zparity };
    Type qubit_type;
};

struct edge_t : graph::base::edge_t {
};

// ID tools that use the above masks. This is only recommended if you
// use the IDs directly.
bool    is_data(uint64_t id) { return id & VERTEX_ID_DATA_FLAG; }
bool    is_xparity(uint64_t id) { return id & VERTEX_ID_XPARITY_FLAG; }
bool    is_zparity(uint64_t id) { return id & VERTEX_ID_ZPARITY_FLAG: }

// Debugging print statement for tanner vertices.
std::string
print_id(uint64_t id) {
    uint64_t q = id & VERTEX_ID_NUMBER_MASK;
    std::string prefix;
    if (is_data(id))            prefix = "d";
    else if (is_xparity(id))    prefix = "x";
    else                        prefix = "z";

    return prefix + std::to_string(id);
}

}   // tanner

#define __TannerGraphParent graph::Graph<tanner::vertex_t, tanner::edge_t>

class TannerGraph : public __TannerGraphParent {
public:
    TannerGraph(void)
        :__TannerGraphParent(), 
        data_qubits(),
        xparity_checks(),
        zparity_checks(),
        x_obs_list(),
        z_obs_list()
    {}

    TannerGraph(const TannerGraph& other)
        :__TannerGraphParent(other),
        data_qubits(other.data_qubits),
        xparity_checks(other.xparity_checks),
        zparity_checks(other.zparity_checks),
        x_obs_list(other.x_obs_list),
        z_obs_list(other.z_obs_list)
    {}

    bool add_vertex(tanner::vertex_t* v) override {
        if (!__TannerGraphParent::add_vertex(v))  return false;
        if (v->qubit_type == tanner::vertex_t::Type::data)      data_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::Type::xparity)   xparity_checks.push_back(v);
        if (v->qubit_type == tanner::vertex_t::Type::zparity)   zparity_checks.push_back(v);
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
        if (v->qubit_type == tanner::vertex_t::Type::xparity)   cat = &xparity_qubits;
        if (v->qubit_type == tanner::vertex_t::Type::zparity)   cat = &zparity_qubits;
        for (auto it = cat->begin(); it != cat->end();) {
            if (*it == v)   it = cat->erase(it);
            else            it++;
        }
        __TannerGraphParent::delete_vertex(v);
    }

    std::vector<tanner::vertex_t*> get_vertices_by_type(tanner::vertex_t::Type type) {
        if (type == tanner::vertex_t::Type::data)           return data_qubits;
        else if (type == tanner::vertex_t::Type::xparity)   return xparity_checks;
        else if (type == tanner::vertex_t::Type::zparity)   return zparity_checks;
    }

    std::vector<tanner::vertex_t*> get_checks() {
        std::vector<tanner::vertex_t*> parity_qubits(xparity_checks);
        for (auto c : zparity_checks)   parity_qubits.push_back(c);
        return parity_qubits;
    }

    std::vector<obs_t> get_obs(bool get_x_obs) {
        if (get_x_obs)  return x_obs_list;
        else            return z_obs_list;
    }

    typedef std::vector<tanner::vertex_t*>  obs_t;
private:
    std::vector<tanner::vertex_t*>  data_qubits;
    std::vector<tanner::vertex_t*>  xparity_checks;
    std::vector<tanner::vertex_t*>  zparity_checks;

    std::vector<obs_t>  x_obs_list;
    std::vector<obs_t>  z_obs_list;
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
