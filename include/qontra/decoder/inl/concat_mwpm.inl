/*
 *  author: Suhas Vittal
 *  date:   29 May 2024
 * */

namespace qontra {

inline std::pair<sptr<gd::vertex_t>, sptr<gd::vertex_t>>
make_ev_pair(sptr<gd::vertex_t> x, sptr<gd::vertex_t> y) {
    if (x > y) return std::make_pair(y,x);
    else        return std::make_pair(x,y);
}

}   // qontra
