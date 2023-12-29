/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef LINPROG_BASE_h
#define LINPROG_BASE_h

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

    bool operator==(const lp_var_t& other) const;
    bool operator<(const lp_var_t<T>& other) const;
};

class lp_expr_t {
public:
    lp_expr_t(double);
    lp_expr_t(lp_var_t);
    lp_expr_t(const lp_expr_t&);
    lp_expr_t(lp_expr_t&&);

    lp_expr_t& operator=(const lp_expr_t&);
    
    // Arithmetic operator overloads:
    lp_expr_t& operator+=(lp_expr_t);
    lp_expr_t& operator-=(lp_expr_t);

    lp_expr_t& operator+=(double);
    lp_expr_t& operator-=(double);
    lp_expr_t& operator*=(double);
    lp_expr_t& operator/=(double);

    double constant;
    std::map<lp_var_t, double> coefs;
private:
    template <class FUNC> void apply_vector_lambda(lp_expr_t, FUNC);
    template <class FUNC> void apply_scalar_lambda(FUNC);
};

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

    lp_expr_t<T> lhs;
    double rhs;
    ConstraintDirection relation;
};

#include "base.inl"

#endif
