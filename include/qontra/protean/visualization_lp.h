/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#ifndef PROTEAN_VISUALIZATION_IMPL_h
#define PROTEAN_VISUALIZATION_IMPL_h

#include "qontra/protean/network.h"
#include "qontra/protean/visualization.h"

#include <vtils/linprog/manager.h>

namespace qontra {
namespace protean {

typedef std::tuple<sptr<net::phys_vertex_t>, char>  lp_label_t;
typedef vtils::CPXLPManager<lp_label_t> LP;

lp_label_t  get_x(sptr<net::phys_vertex_t>);
lp_label_t  get_y(sptr<net::phys_vertex_t>);

vtils::lp_expr_t    get_distance_metric(vtils::lp_var_t, vtils::lp_var_t, vtils::lp_var_t, vtils::lp_var_t);

// LP formulation functions:
void    lp_add_variables(LP&, PhysicalNetwork&);

// Constraint: each vertex should be at least D units away.
//
// Formally:
//      (xi - xj)^2 + (yi - yj)^2 >= D**2
//  --> xi**2 - 2*xi*xj + xj**2 + ... >= D**2
void    lp_add_minimum_distance_constraints(LP&, PhysicalNetwork&, placement_config_t);

// Create constraints to identify crossing edges.
// Let AB and CD be line segments. We say AB and CD are x-crossing if
//      (xA - xC)(xA - xD)(xB - xC)(xB - xD) > 0
// and similarly for y-crossing. AB and CD are intersecting if they are both
// x-crossing and y-crossing.
//
// We can check for x-crossing and y-crossing via linear programming. Let
// bX_1 = 1 if xA - xC > 0 (and etc. for bX_2 ... bX_4). Do the same for bY_i.
//
// Let BX = 1 if the sum of bX_i's is even.
//      --> sum(bX_i) = 2y + BX - 1 where 0 <= y <= 2 is an integer.
// AB and CD are intersecting if BX + BY = 2.
//      --> 2*ind + 1 >= BX + BY where ind = 1 if AB and CD are intersecting.
vtils::lp_expr_t    lp_add_crossing_edges_objective(LP&, PhysicalNetwork&, placement_config_t);
vtils::lp_expr_t    lp_add_edge_distance_objective(LP&, PhysicalNetwork&);

}   // protean
}   // qontra

#include "visualization_lp.inl"

#endif  // PROTEAN_VISUALIZATION_IMPL_h
