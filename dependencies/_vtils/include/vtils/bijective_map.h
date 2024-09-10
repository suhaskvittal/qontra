/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#ifndef VTILS_BIJECTIVE_MAP_h
#define VTILS_BIJECTIVE_MAP_h

#include <unordered_map>
#include <type_traits>

namespace vtils {

template <class T, class U, template<class,class> class MAP_TYPE> 
MAP_TYPE<U, T> get_reverse_map(const MAP_TYPE<T, U>&);

// This is for the common situation where you just need to label something
// with another thing. This class enforces bijectivity.

#define TYPE_ASSERT(fn) static_assert(!std::is_same<T, U>::value, "function #fn is undefined if T = U.")

template < class T, class U, template<class,class> class MAP_TYPE=std::unordered_map >
class BijectiveMap {
public:
    BijectiveMap()
        :forward_map(),
        reverse_map()
    {}

    BijectiveMap(const MAP_TYPE<T, U>& orig)
        :forward_map(orig),
        reverse_map(get_reverse_map(orig))
    {}

    BijectiveMap(const BijectiveMap& other)
        :forward_map(other.forward_map),
        reverse_map(other.reverse_map)
    {}

    BijectiveMap(BijectiveMap&& other)
        :forward_map(std::move(other.forward_map)),
        reverse_map(std::move(other.reverse_map))
    {}

    inline BijectiveMap& operator=(const BijectiveMap& other);

    BijectiveMap reverse(void);

    inline void clear(void) { forward_map.clear(); reverse_map.clear(); }

    // Pretty much all the functions below are wrappers for MAP_TYPE.
    
    // Returns false if change is not bijective.
    inline bool put(T key, U value) { TYPE_ASSERT("put"); return f_put(key, value); }
    inline bool put(U key, T value) { TYPE_ASSERT("put"); return r_put(key, value); }

    inline U at(T key) const { TYPE_ASSERT("at"); return f_at(key); }
    inline T at(U key) const { TYPE_ASSERT("at"); return r_at(key); }

    inline void erase(T key) { TYPE_ASSERT("erase"); f_erase(key); }
    inline void erase(U key) { TYPE_ASSERT("erase"); r_erase(key); }

    inline size_t count(T key) const { TYPE_ASSERT("count"); return f_count(key); }
    inline size_t count(U key) const { TYPE_ASSERT("count"); return r_count(key); }

    inline void swap(T x, T y) { TYPE_ASSERT("swap"); f_swap(x, y); }
    inline void swap(U x, U y) { TYPE_ASSERT("swap"); r_swap(x, y); }

    // Specific functions for interacting with the forward and reverse maps:
    // 
    // This is useful when T = U as the above functions are deleted.

    bool f_put(T, U);
    bool r_put(U, T);

    void f_erase(T);
    void r_erase(U);

    inline U f_at(T key) const { return forward_map.at(key); }
    inline T r_at(U key) const { return reverse_map.at(key); }

    inline size_t f_count(T key) const { return forward_map.count(key); }
    inline size_t r_count(U key) const { return reverse_map.count(key); }

    void f_swap(T, T);  // Swaps the values for the two keys.
    void r_swap(U, U);

    // Other functions:
    inline MAP_TYPE<T, U>::iterator begin() { return forward_map.begin(); }
    inline MAP_TYPE<T, U>::const_iterator cbegin() const { return forward_map.cbegin(); }
    inline MAP_TYPE<T, U>::iterator end() { return forward_map.end(); }
    inline MAP_TYPE<T, U>::const_iterator cend() const { return forward_map.cend(); }

    inline size_t size(void) { return forward_map.size(); }

    bool assert_bijectivity = true;
private:
    static BijectiveMap from_maps(const MAP_TYPE<T, U>& fwd, const MAP_TYPE<U, T>& rev) {
        BijectiveMap bmap;
        bmap.forward_map = fwd;
        bmap.reverse_map = rev;
        return bmap;
    }

    MAP_TYPE<T, U> forward_map;
    MAP_TYPE<U, T> reverse_map;
};

}   // vtils

#include "inl/bijective_map.inl"

#endif  // DEFS_BIJECTIVE_MAP_h
