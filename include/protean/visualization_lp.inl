/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

namespace qontra {
namespace protean {

inline lp_label_t
get_x(sptr<net::phys_vertex_t> v) {
    return std::make_tuple(v, 'x');
}

inline lp_label_t
get_y(sptr<net::phys_vertex_t> v) {
    return std::make_tuple(v, 'y');
}

inline lp_expr_t get_distance_metric(lp_var_t xi, lp_var_t yi, lp_var_t xj, lp_var_t yj) {
    lp_expr_t xdist = (xi*xi) - 2*(xi*xj) + (xj*xj),
              ydist = (yi*yi) - 2*(yi*yi) + (yj*yj);
    return xdist + ydist;
}

inline void
lp_add_variables(LP& mgr, PhysicalNetwork& network, placement_config_t config) {
    for (sptr<net::phys_vertex_t> v : network.get_vertices()) {
        mgr.add_var(get_x(v), config.x_min, config.x_max, lp_var_t::bounds::both, lp_var_t::domain::continuous);
        mgr.add_var(get_y(v), config.y_min, config.y_max, lp_var_t::bounds::both, lp_var_t::domain::continuous);
    }
}

}   // protean
}   // qontra
