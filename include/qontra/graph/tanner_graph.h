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

class TannerGraph;

namespace io {

// Function for io callback.
//
// In a tanner graph description file, each line should be of the form:
//          <X/Z><check-id>,<data-qubit-1>,<data-qubit-2>,...

void update_tanner_graph(TannerGraph&, std::string); // Callback for io function.

}   // io

class TannerGraph : public Graph<tanner::vertex_t, tanner::edge_t> {
public:
    typedef std::vector<sptr<tanner::vertex_t>>  obs_t;

    TannerGraph(void) = default;
    TannerGraph(TannerGraph&& other) = default;

    bool    add_vertex(sptr<tanner::vertex_t>) override;
    bool    add_edge(sptr<tanner::edge_t>) override;
    void    delete_vertex(sptr<tanner::vertex_t>) override;
    
    std::vector<sptr<tanner::vertex_t>> get_vertices_by_type(tanner::vertex_t::type) const;
    std::vector<sptr<tanner::vertex_t>> get_checks(void) const;

    // Sets the map reference to the color map. Returns the max color used.
    int compute_check_color_map(std::map<sptr<tanner::vertex_t>, int>&) const;
    int compute_code_distance(bool for_x) const;

    std::vector<obs_t> get_obs(bool get_x_obs) const;
private:
    std::vector<sptr<tanner::vertex_t>>& get_vertices_by_type_(tanner::vertex_t::type);

    // This is a helper function for compute_check_color_map.
    int update_check_color_map(std::map<sptr<tanner::vertex_t>, int>&, bool use_x_checks) const;

    std::vector<sptr<tanner::vertex_t>>  data_qubits;
    std::vector<sptr<tanner::vertex_t>>  xparity_checks;
    std::vector<sptr<tanner::vertex_t>>  zparity_checks;

    std::vector<obs_t>  x_obs_list;
    std::vector<obs_t>  z_obs_list;

    friend void io::update_tanner_graph(TannerGraph&, std::string);
};

// Specialization of print_v to tanner::vertex_t
template <>
std::string print_v(sptr<tanner::vertex_t>);

}   // graph
}   // qontra

#include "tanner_graph.inl"

#endif  // QONTRA_TANNER_GRAPH_h
