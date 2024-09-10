/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

namespace vtils {

inline bool
lp_var_t::operator==(const lp_var_t& other) const {
    return column == other.column;
}

inline bool
lp_var_t::operator<(const lp_var_t& other) const {
    return column < other.column;
}

inline bool
lp_quad_var_t::operator==(const lp_quad_var_t& other) const {
    return v1 == other.v1 && v2 == other.v2;
}

inline bool
lp_quad_var_t::operator<(const lp_quad_var_t& other) const {
    return v1 < other.v1 || (v1 == other.v1 && v2 < other.v2);
}

}   // vtils
