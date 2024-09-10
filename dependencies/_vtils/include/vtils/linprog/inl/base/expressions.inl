/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

namespace vtils {

inline bool
lp_expr_t::is_quadratic() {
    return !q_coefs.empty();
}

inline lp_expr_t&
lp_expr_t::operator=(const lp_expr_t& other) {
    l_coefs = other.l_coefs;
    q_coefs = other.q_coefs;
    constant = other.constant;
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator+=(lp_expr_t e) {
    apply_vector_lambda(e, [] (double x, double y) { return x + y; });
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator-=(lp_expr_t e) {
    apply_vector_lambda(e, [] (double x, double y) { return x - y; });
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator+=(double x) {
    constant += x;
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator-=(double x) {
    constant -=x;
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator*=(double x) {
    apply_scalar_lambda([&] (double a) { return a * x; });
    return *this;
}

inline lp_expr_t&
lp_expr_t::operator/=(double x) {
    apply_scalar_lambda([&] (double a) { return a / x; });
    return *this;
}

template <class FUNC> inline void
lp_expr_t::apply_vector_lambda(lp_expr_t e, FUNC f) {
    constant = f(constant, e.constant);
    for (auto& kv : e.l_coefs) {
        l_coefs[kv.first] = f(l_coefs[kv.first], kv.second);
    }
    for (auto& kv : e.q_coefs) {
        q_coefs[kv.first] = f(q_coefs[kv.first], kv.second);
    }
}

template <class FUNC> inline void
lp_expr_t::apply_scalar_lambda(FUNC f) {
    constant = f(constant);
    for (auto& kv : l_coefs) kv.second = f(kv.second);
    for (auto& kv : q_coefs) kv.second = f(kv.second);
}

}   // vtils

