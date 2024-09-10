/*
 *  author: Suhas Vittal
 *  date:   29 Decemeber 2023
 * */

namespace vtils {

inline lp_constr_t&
lp_constr_t::operator=(const lp_constr_t& other) {
    lhs = other.lhs;
    rhs = other.rhs;
    relation = other.relation;
    return *this;
}

inline bool
lp_constr_t::is_quadratic() {
    return lhs.is_quadratic();
}

}   // vtils
