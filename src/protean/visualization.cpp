/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "protean/visualization.h"

#include <linprog/manager.h>

#include <graphviz/gvc.h>

namespace qontra {
namespace protean {

std::vector<Plane> place_network(PhysicalNetwork& network, placement_config_t config) {
    CPXLPManager<sptr<net::phys_vertex_t>> mgr;
    // Create an LP variable for each qubit in the network.
    std::map<sptr<net::phys_vertex_t>, lp_var_t> vertex_to_x_var;
    std::map<sptr<net::phys_vertex_t>, lp_var_t> vertex_to_y_var;

    for (sptr<net::phys_vertex_t> v : network.get_vertices()) {
        vertex_to_x_var[v] = 
            mgr.add_var(v, config.x_min, config.x_max, lp_var_t::bounds::both, lp_var_t::domain::continuous);
        vertex_to_y_var[v] = 
            mgr.add_var(v, config.y_min, config.y_max, lp_var_t::bounds::both, lp_var_t::domain::continuous);
    }
    // Constraint 1: each vertex should be at least D units away.
    //
    // Formally:
    //      (xi - xj)^2 + (yi - yj)^2 >= D
    //  --> 
}

}   // protean
}   // qontra
