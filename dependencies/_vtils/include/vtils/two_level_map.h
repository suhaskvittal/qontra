/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef VTILS_TWO_LEVEL_MAP_h
#define VTILS_TWO_LEVEL_MAP_h

#include <unordered_map>

namespace vtils {

template < class T, class U, class V, template<class,class> class MAP_TYPE=std::unordered_map >
using TwoLevelMap = MAP_TYPE<T, MAP_TYPE<U, V>>;

template <class T, class U, class V> inline void
tlm_put(TwoLevelMap<T, U, V>& m, T x, U y, V z) {
    m[x][y] = z;
}

template <class T, class U, class V> inline V
tlm_at(const TwoLevelMap<T, U, V>& m, T x, U y) {
    return m.at(x).at(y);
}

}   // vtils

#endif  // VTILS_TWO_LEVEL_MAP_h
