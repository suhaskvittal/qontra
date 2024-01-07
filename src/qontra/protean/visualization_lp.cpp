/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "qontra/protean/visualization_lp.h"

using namespace vtils;

namespace qontra {
namespace protean {

inline fp_t sqr(fp_t x) {
    return x*x;
}

void
lp_add_minimum_distance_constraints(LP& mgr, PhysicalNetwork& network, placement_config_t config) {
    auto vertices = network.get_vertices();
    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<net::phys_vertex_t> v = vertices[i];
        lp_var_t xi = mgr.get_var(get_x(v)),
                 yi = mgr.get_var(get_y(v));
        for (size_t j = i+1; j < vertices.size(); j++) {
            sptr<net::phys_vertex_t> w = vertices[j];
            lp_var_t xj = mgr.get_var(get_x(w)),
                     yj = mgr.get_var(get_y(w));
            // Make constraint.
            lp_constr_t con(get_distance_metric(xi, yi, xj, yj), 
                            sqr(config.min_distance_between_qubits),
                            lp_constr_t::direction::ge);
            mgr.add_constraint(con);
        }
    }
}

lp_expr_t
lp_add_crossing_edges_objective(LP& mgr, PhysicalNetwork& network, placement_config_t config) {
    const fp_t M = config.x_max - config.x_min;
    const lp_var_t::domain svardom = config.edge_crossing_relax_variables ? 
                                        lp_var_t::domain::continuous : lp_var_t::domain::integer;

    lp_expr_t crossing_edge_sum;
    auto edges = network.get_edges();

    size_t indicator_count = 0;
    for (size_t i = 0; i < edges.size(); i++) {
        sptr<net::phys_edge_t> e1 = edges[i];
        if (e1->is_out_of_plane()) continue;

        sptr<net::phys_vertex_t> v11 = std::reinterpret_pointer_cast<net::phys_vertex_t>(e1->src),
                                 v12 = std::reinterpret_pointer_cast<net::phys_vertex_t>(e1->dst);
        lp_var_t x11 = mgr.get_var(get_x(v11)),
                 y11 = mgr.get_var(get_y(v11)),
                 x12 = mgr.get_var(get_x(v12)),
                 y12 = mgr.get_var(get_y(v12));
        for (size_t j = i+1; j < edges.size(); j++) {
            sptr<net::phys_edge_t> e2 = edges[j];
            if (e2->is_out_of_plane()) continue;

            sptr<net::phys_vertex_t> v21 = std::reinterpret_pointer_cast<net::phys_vertex_t>(e2->src),
                                     v22 = std::reinterpret_pointer_cast<net::phys_vertex_t>(e2->dst);
            // Skip if any endpoints are shared.
            if (v11 == v21 || v11 == v22 || v12 == v21 || v12 == v22) continue;

            lp_var_t x21 = mgr.get_var(get_x(v21)),
                     y21 = mgr.get_var(get_y(v21)),
                     x22 = mgr.get_var(get_x(v22)),
                     y22 = mgr.get_var(get_y(v22));
            // Make the slack variables to detect crossing edges.
            std::array<lp_var_t, 4> xi_var_array{x11, x11, x12, x12},
                                    xj_var_array{x21, x22, x21, x22},
                                    yi_var_array{y11, y11, y12, y12},
                                    yj_var_array{y21, y22, y21, y22};
            lp_expr_t bx_sum, by_sum;
            for (size_t i = 0; i < 4; i++) {
                lp_var_t xi = xi_var_array[i],
                         xj = xj_var_array[i],
                         yi = yi_var_array[i],
                         yj = yj_var_array[i];

                lp_var_t bx = mgr.add_slack_var(0.0, 1.0, lp_var_t::bounds::both, svardom),
                         by = mgr.add_slack_var(0.0, 1.0, lp_var_t::bounds::both, svardom);
                // x1c constrains bx = 1 if xi - xj > 0
                // x2c constrains bx = 0 if xi - xj < 0.
                // similarly for y1c and y2c.
                lp_constr_t x1c(xi - xj - M*bx, 0, lp_constr_t::direction::le),
                            x2c(xi - xj + M*bx, M, lp_constr_t::direction::le),
                            y1c(yi - yj - M*by, 0, lp_constr_t::direction::le),
                            y2c(yi - yj + M*by, M, lp_constr_t::direction::le);
                mgr.add_constraint(x1c);
                mgr.add_constraint(x2c);
                mgr.add_constraint(y1c);
                mgr.add_constraint(y2c);
                bx_sum += bx;
                by_sum += by;
            }
            // Now, make BX (bxx) and BY (byy). xs and ys are other slack variables.
            lp_var_t bxx = mgr.add_slack_var(0.0, 1.0, lp_var_t::bounds::both, svardom),
                     byy = mgr.add_slack_var(0.0, 1.0, lp_var_t::bounds::both, svardom),
                     xs = mgr.add_slack_var(0.0, 2.0, lp_var_t::bounds::both, lp_var_t::domain::integer),
                     ys = mgr.add_slack_var(0.0, 2.0, lp_var_t::bounds::both, lp_var_t::domain::integer);
            // The below constraints bxx (byy) to be 1 if bx_sum (by_sum) is even.
            lp_constr_t bxxc(2*xs + bxx - 1 - bx_sum, 0, lp_constr_t::direction::eq);
            lp_constr_t byyc(2*ys + byy - 1 - bx_sum, 0, lp_constr_t::direction::eq);
            mgr.add_constraint(bxxc);
            mgr.add_constraint(byyc);
            // Finally, make the constraint for the indicator variable.
            lp_var_t ind = mgr.add_slack_var(0.0, 1.0, lp_var_t::bounds::both, svardom);
            lp_constr_t indc(2*ind + 1 - bxx - byy, 0, lp_constr_t::direction::ge);
            mgr.add_constraint(indc);

            crossing_edge_sum += ind;

            if (++indicator_count == config.edge_crossing_max_indicators) goto crossing_edge_exit;
        }
    }
crossing_edge_exit:
    return crossing_edge_sum;
}

lp_expr_t
lp_add_edge_distance_objective(LP& mgr, PhysicalNetwork& network) {
    lp_expr_t edge_length_sum;
    for (sptr<net::phys_vertex_t> v : network.get_vertices()) {
        lp_var_t xi = mgr.get_var(get_x(v)),
                 yi = mgr.get_var(get_y(v));
        auto neighbors = network.get_neighbors(v);
        for (size_t i = 0; i < neighbors.size(); i++) {
            sptr<net::phys_vertex_t> w = neighbors[i];
            sptr<net::phys_edge_t> e = network.get_edge(v, w);
            if (e->is_out_of_plane()) continue;
            // Otherwise, we can continue.
            lp_var_t xj = mgr.get_var(get_x(w)),
                     yj = mgr.get_var(get_y(w));
            edge_length_sum += get_distance_metric(xi, yi, xj, yj);
        }
    }
    return edge_length_sum;
}

}   // protean
}   // qontra
