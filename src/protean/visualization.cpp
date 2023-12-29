/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "protean/visualization_impl.h"

#include <graphviz/gvc.h>

#include <iostream>

namespace qontra {
namespace protean {

inline fp_t sqr(fp_t x) {
    return x*x;
}

Plane
place_network(PhysicalNetwork& network, placement_config_t config) {
    // Build the LP.
    LP mgr;
    lp_add_variables(mgr, network);
    lp_add_minimum_distance_constraints(mgr, network, config);

    lp_expr_t objective;
    if (!config.edge_crossing_skip) {
        objective += 
            config.edge_crossing_objective_coef * lp_add_crossing_edges_objective(mgr, network, config);
    }
    objective += config.edge_length_objective_coef * lp_add_edge_distance_objective(mgr, network);
    mgr.build(objective, false);
    // Now, we can get the result and build the plane.
    fp_t obj;
    int solstat;
    if (mgr.solve(&obj, &solstat)) {
        std::cerr << "place_network: program is infeasible" << std::endl;
        return Plane();
    }
    // Form layout from the LP results.
    Plane layout;
    for (sptr<net::phys_vertex_t> v : network.get_vertices()) {
        coord_t<2> p{ mgr.get_value(get_x(v)), mgr.get_value(get_y(v)) };
        layout.put(v, p);
    }
    return layout;
}

}   // protean
}   // qontra
