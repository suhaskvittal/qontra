/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

lp_expr_t::lp_expr_t(double c=0.0)
    :coefs(),
    constant(c)
{}

lp_expr_t::lp_expr_t(lp_var_t v)
    :coefs(),
    constant(0)
{
    coefs[v] = 1.0;
}

lp_expr_t::lp_expr_t(const lp_expr_t& other)
    :coefs(other.coefs),
    constant(other.constant)
{}

lp_expr_t::lp_expr_t(lp_expr_t&& other)
    :coefs(std::move(other.coefs)),
    constant(other.constant)
{}

inline lp_expr_t&
lp_expr_t::operator=(const lp_expr_t& other) {
    coefs = other.coefs;
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
    for (auto& kv : e.coefs) {
        coefs[kv.first] = f(coefs[kv.first], kv.second);
    }
}

template <class FUNC> inline void
lp_expr_t::apply_scalar_lambda(FUNC f) {
    constant = f(constant);
    for (auto& kv : coefs) {
        kv.second = f(kv.second);
    }
}

inline lp_expr_t operator-(lp_expr_t e) {
    e.constant = -e.constant;
    for (auto& kv : e.coefs) kv.second = -kv.second;
    return e;
}

inline lp_expr_t operator+(lp_expr_t a, lp_expr_t b) { a += b; return a; }
inline lp_expr_t operator-(lp_expr_t a, lp_expr_t b) { a -= b; return a; }

inline lp_expr_t operator+(lp_expr_t a, double x) { a += x; return a; }
inline lp_expr_t operator-(lp_expr_t a, double x) { a -= x; return a; }
inline lp_expr_t operator*(lp_expr_t a, double x) { a *= x; return a; }
inline lp_expr_t operator/(lp_expr_t a, double x) { a /= x; return a; }

inline lp_expr_t operator+(double x, lp_expr_t a) { return a + x; }
inline lp_expr_t operator-(double x, lp_expr_t a) { return (-a) + x; }
inline lp_expr_t operator*(double x, lp_expr_t a) { return a * x; }

inline bool
lp_var_t::operator==(const lp_var_t& other) const {
    return column == other.column;
}

inline bool
lp_var_t::operator<(const lp_var_t& other) const {
    return column < other.column;
}

lp_constr_t::lp_constr_t(lp_expr_t _lhs, lp_expr_t _rhs, lp_constr_t::direction r)
    :lhs(_lhs),
    rhs(_rhs),
    relation(r)
{
    // Move rhs to the lhs.
    for (auto kv : rhs.coefs) {
        lhs.coefs[kv.first] -= kv.second;
    }
    rhs = static_cast<lp_expr_t>(rhs.constant - lhs.constant);
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

