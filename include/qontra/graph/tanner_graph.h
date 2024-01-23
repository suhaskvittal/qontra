/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#ifndef QONTRA_TANNER_GRAPH_h
#define QONTRA_TANNER_GRAPH_h

#include "qontra/defs.h"
#include "qontra/graph.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace qontra {
namespace graph {

namespace tanner {

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
    enum class type { data, xparity, zparity };
    type qubit_type;
};

struct edge_t : graph::base::edge_t {
};

// ID tools that use the above masks. This is only recommended if you
// use the IDs directly.
inline bool is_data(uint64_t id) { return id & VERTEX_ID_DATA_FLAG; }
inline bool is_xparity(uint64_t id) { return id & VERTEX_ID_XPARITY_FLAG; }
inline bool is_zparity(uint64_t id) { return id & VERTEX_ID_ZPARITY_FLAG; }

// Debugging print statement for tanner vertices.
inline std::string
print_id(uint64_t id) {
    uint64_t q = id & VERTEX_ID_NUMBER_MASK;
    std::string prefix;
    if (is_data(id))            prefix = "d";
    else if (is_xparity(id))    prefix = "x";
    else                        prefix = "z";

    return prefix + std::to_string(q);
}

}   // tanner

#define __TannerGraphParent graph::Graph<tanner::vertex_t, tanner::edge_t>

class TannerGraph;

namespace io {

// Function for io callback.
//
// In a tanner graph description file, each line should be of the form:
//          <X/Z><check-id>,<data-qubit-1>,<data-qubit-2>,...

void update_tanner_graph(TannerGraph&, std::string); // Callback for io function.

}   // io

class TannerGraph : public __TannerGraphParent {
public:
    typedef std::vector<sptr<tanner::vertex_t>>  obs_t;

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

    bool add_vertex(sptr<tanner::vertex_t> v) override {
        if (!__TannerGraphParent::add_vertex(v))  return false;
        if (v->qubit_type == tanner::vertex_t::type::data)      data_qubits.push_back(v);
        if (v->qubit_type == tanner::vertex_t::type::xparity)   xparity_checks.push_back(v);
        if (v->qubit_type == tanner::vertex_t::type::zparity)   zparity_checks.push_back(v);
        return true;
    }

    bool add_edge(sptr<tanner::edge_t> e) override {
        auto v = std::reinterpret_pointer_cast<tanner::vertex_t>(e->src);
        auto w = std::reinterpret_pointer_cast<tanner::vertex_t>(e->dst);
        // Make sure the edge preserves the bipartite property.
        bool src_is_parity = v->qubit_type != tanner::vertex_t::type::data;
        bool dst_is_parity = w->qubit_type != tanner::vertex_t::type::data;
        if (src_is_parity && dst_is_parity) return false;
        return __TannerGraphParent::add_edge(e);
    }

    void delete_vertex(sptr<tanner::vertex_t> v) override {
        std::vector<sptr<tanner::vertex_t>>* cat;
        if (v->qubit_type == tanner::vertex_t::type::data)      cat = &data_qubits;
        if (v->qubit_type == tanner::vertex_t::type::xparity)   cat = &xparity_checks;
        if (v->qubit_type == tanner::vertex_t::type::zparity)   cat = &zparity_checks;
        for (auto it = cat->begin(); it != cat->end();) {
            if (*it == v)   it = cat->erase(it);
            else            it++;
        }
        __TannerGraphParent::delete_vertex(v);
    }

    std::vector<sptr<tanner::vertex_t>> get_vertices_by_type(tanner::vertex_t::type type) {
        if (type == tanner::vertex_t::type::data)           return data_qubits;
        else if (type == tanner::vertex_t::type::xparity)   return xparity_checks;
        else                                                return zparity_checks;
    }

    std::vector<sptr<tanner::vertex_t>> get_checks() {
        std::vector<sptr<tanner::vertex_t>> parity_qubits(xparity_checks);
        for (auto c : zparity_checks)   parity_qubits.push_back(c);
        return parity_qubits;
    }

    std::vector<obs_t> get_obs(bool get_x_obs) {
        if (get_x_obs)  return x_obs_list;
        else            return z_obs_list;
    }
private:
    std::vector<sptr<tanner::vertex_t>>  data_qubits;
    std::vector<sptr<tanner::vertex_t>>  xparity_checks;
    std::vector<sptr<tanner::vertex_t>>  zparity_checks;

    std::vector<obs_t>  x_obs_list;
    std::vector<obs_t>  z_obs_list;

    friend void io::update_tanner_graph(TannerGraph&, std::string);
};

// Specialization of print_v to tanner::vertex_t
template <> inline std::string
print_v(sptr<tanner::vertex_t> v) {
    std::string type_prefix;
    if (v->qubit_type == tanner::vertex_t::type::data) {
        type_prefix = "d";
    } else if (v->qubit_type == tanner::vertex_t::type::xparity) {
        type_prefix = "x";
    } else {
        type_prefix = "z";
    }
    return type_prefix + std::to_string(v->id & tanner::VERTEX_ID_NUMBER_MASK);
}

}   // graph
}   // qontra

#endif  // QONTRA_TANNER_GRAPH_h
