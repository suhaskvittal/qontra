/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef DEFS_TWO_LEVEL_MAP_h
#define DEFS_TWO_LEVEL_MAP_h

#include <map>

template <class T, class U, class V> using TwoLevelMap = std::map<T, std::map<U, V>>;

template <class T, class U, class V> inline void
tlm_put(TwoLevelMap<T, U, V>& m, T x, U y, V z) {
    if (!m.count(x))    m[x] = std::map<U, V>();
    m[x][y] = z;
}

#endif  // DEFS_TWO_LEVEL_MAP_h
