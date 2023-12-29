/*
 *  author: Suhas Vittal
 *  date:   29 Decemeber 2023
 * */

lp_constr_t::lp_constr_t(lp_expr_t _lhs, lp_expr_t _rhs, lp_constr_t::direction r)
    :lhs(_lhs),
    rhs(_rhs),
    relation(r)
{
    // Move rhs to the lhs.
    lhs -= rhs;
    rhs = static_cast<lp_expr_t>(-lhs.constant);
    lhs.constant = 0;
}

lp_constr_t::lp_constr_t(lp_expr_t _lhs, double _rhs, lp_constr_t::direction r)
    :lhs(_lhs),
    rhs(_rhs),
    relation(r)
{
    rhs -= lhs.constant;
    lhs.constant = 0;
}

inline bool
lp_constr_t::is_quadratic() {
    return lhs.q_coefs().size() > 0;
}

