/*
 *  author: Suhas Vittal
 *  date:   2 January 2024
 * */

#include "vtils/linprog/base.h"

#include <utility>

namespace vtils {

lp_expr_t::lp_expr_t(double c)
    :l_coefs(),
    q_coefs(),
    constant(c)
{}

lp_expr_t::lp_expr_t(lp_var_t v)
    :l_coefs(),
    q_coefs(),
    constant(0)
{
    l_coefs[v] = 1.0;
}

lp_expr_t::lp_expr_t(lp_quad_var_t v)
    :l_coefs(),
    q_coefs(),
    constant(0)
{
    q_coefs[v] = 1.0;
}

lp_expr_t::lp_expr_t(const lp_expr_t& other)
    :l_coefs(other.l_coefs),
    q_coefs(other.q_coefs),
    constant(other.constant)
{}

lp_expr_t::lp_expr_t(lp_expr_t&& other)
    :l_coefs(std::move(other.l_coefs)),
    q_coefs(std::move(other.q_coefs)),
    constant(other.constant)
{}

lp_constr_t::lp_constr_t(lp_expr_t _lhs, lp_expr_t _rhs, lp_constr_t::direction r)
    :lhs(),
    rhs(),
    relation(r)
{
    // Move rhs to the lhs.
    _lhs -= _rhs;
    lhs = _lhs;
    rhs = -_lhs.constant;
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

lp_constr_t::lp_constr_t(const lp_constr_t& other)
    :lhs(other.lhs),
    rhs(other.rhs),
    relation(other.relation)
{}

lp_constr_t::lp_constr_t(lp_constr_t&& other)
    :lhs(std::move(other.lhs)),
    rhs(std::move(other.rhs)),
    relation(other.relation)
{}

}   // vtils
