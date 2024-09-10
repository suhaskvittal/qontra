/*
 *  author: Suhas Vittal
 *  
 *  Hashing extensions for commonly used types.
 * */

#ifndef VTILS_HASHEXT_h
#define VTILS_HASHEXT_h

#include <bit>
#include <utility>
// Need a header where hash is defined:
#include <unordered_set>

namespace std {

template <class T, class U>
struct hash< std::pair<T,U> > {
    using argument_type=pair<T,U>;
    using result_type = std::size_t;
        
    result_type operator()(const argument_type& p) const {
        size_t h1 = hash<T>{}(p.first),
               h2 = hash<U>{}(p.second);
        return rotl(h1, 8) ^ h2;
    }
};

}   // std

#endif  // VTILS_HASHEXT_h
