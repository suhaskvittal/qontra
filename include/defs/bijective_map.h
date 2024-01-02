/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#ifndef DEFS_BIJECTIVE_MAP_h
#define DEFS_BIJECTIVE_MAP_h

#include <map>
#include <string>

template <class T, class U> 
std::map<U, T> get_reverse_map(const std::map<T, U>&);

// This is for the common situation where you just need to label something
// with another thing. This class enforces bijectivity.

#define TYPE_ASSERT(fn) static_assert(!std::is_same<T, U>::value, "function #fn is undefined if T = U.")

template <class T, class U>
class BijectiveMap {
public:
    BijectiveMap()
        :forward_map(),
        reverse_map()
    {}

    BijectiveMap(const std::map<T, U>& orig)
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

    // Pretty much all the functions below are wrappers for std::map.
    
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
    inline std::map<T, U>::iterator begin() { return forward_map.begin(); }
    inline std::map<T, U>::const_iterator cbegin() const { return forward_map.cbegin(); }
    inline std::map<T, U>::iterator end() { return forward_map.end(); }
    inline std::map<T, U>::const_iterator cend() const { return forward_map.cend(); }

    inline size_t size(void) { return forward_map.size(); }

    bool assert_bijectivity = true;
private:
    static BijectiveMap from_maps(const std::map<T, U>& fwd, const std::map<U, T>& rev) {
        BijectiveMap bmap;
        bmap.forward_map = fwd;
        bmap.reverse_map = rev;
        return bmap;
    }

    std::map<T, U> forward_map;
    std::map<U, T> reverse_map;
};

#include "bijective_map.inl"

#endif  // DEFS_BIJECTIVE_MAP_h
