/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef LINPROG_BASE_h
#define LINPROG_BASE_h

#include <map>

#include <stddef.h>

struct lp_var_t {
    // Variable id:
    size_t  column;
    // Bounds:
    double  lwr;
    double  upp;

    enum class bounds { none, lower, upper, both, fixed };
    enum class domain { continuous, integer, binary };

    bounds bounds_type;
    domain var_type;

    bool operator==(const lp_var_t&) const;
    bool operator<(const lp_var_t&) const;
};

struct lp_quad_var_t {
    lp_var_t    v1;
    lp_var_t    v2;

    bool operator==(const lp_quad_var_t&) const;
    bool operator<(const lp_quad_var_t&) const;
};

struct lp_expr_t {
    lp_expr_t(double=0.0);
    lp_expr_t(lp_var_t);
    lp_expr_t(lp_quad_var_t);
    lp_expr_t(const lp_expr_t&);
    lp_expr_t(lp_expr_t&&);

    bool is_quadratic(void);

    lp_expr_t& operator=(const lp_expr_t&);
    
    // Arithmetic operator overloads:
    lp_expr_t& operator+=(lp_expr_t);
    lp_expr_t& operator-=(lp_expr_t);

    lp_expr_t& operator+=(double);
    lp_expr_t& operator-=(double);
    lp_expr_t& operator*=(double);
    lp_expr_t& operator/=(double);

    double constant;
    std::map<lp_var_t, double>      l_coefs;
    std::map<lp_quad_var_t, double> q_coefs;

    template <class FUNC> void apply_vector_lambda(lp_expr_t, FUNC);
    template <class FUNC> void apply_scalar_lambda(FUNC);
};

// All operator overloads:
//
// The product operator for variables is special:
lp_quad_var_t operator*(lp_var_t, lp_var_t);

lp_expr_t operator-(lp_expr_t);

lp_expr_t operator+(lp_expr_t, lp_expr_t);
lp_expr_t operator-(lp_expr_t, lp_expr_t);

lp_expr_t operator+(lp_expr_t, double);
lp_expr_t operator-(lp_expr_t, double);
lp_expr_t operator*(lp_expr_t, double);
lp_expr_t operator/(lp_expr_t, double);

lp_expr_t operator+(double, lp_expr_t);
lp_expr_t operator-(double, lp_expr_t);
lp_expr_t operator*(double, lp_expr_t);

// NOTE: left division is not defined as this is nonlinear.

struct lp_constr_t {
    enum class direction { ge, le, eq, neq };

    lp_constr_t(lp_expr_t, lp_expr_t, direction);
    lp_constr_t(lp_expr_t, double, direction);

    lp_constr_t(const lp_constr_t&);
    lp_constr_t(lp_constr_t&&);

    lp_constr_t& operator=(const lp_constr_t&);

    bool is_quadratic(void);

    lp_expr_t   lhs;
    double      rhs;
    direction   relation;
};

#include "variables.inl"
#include "expressions.inl"
#include "constraints.inl"
#include "operations.inl"

#endif  // LINPROG_BASE_h
