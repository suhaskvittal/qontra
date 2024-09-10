/*
 *  author: Suhas Vittal
 *  date:   29 Decemeber 2023
 * */

namespace vtils {

inline lp_quad_var_t operator*(lp_var_t v, lp_var_t w) {
    lp_quad_var_t q;
    q.v1 = v;
    q.v2 = w;
    return q;
} 

inline lp_expr_t operator-(lp_expr_t e) {
    e.apply_scalar_lambda([] (double x) { return -x; });
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

}   // vtils
