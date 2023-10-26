/* 
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DEFS_h
#define DEFS_h

#include <stim.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <map>
#include <vector>

#include <stdint.h>

typedef double fp_t;

#ifdef __APPLE__
typedef uint32_t uint;
#endif

typedef uint64_t                addr_t;
typedef std::vector<uint64_t>   vlw_t;
typedef stim::simd_bits         syndrome_t;

const fp_t KB = 1024.0;
const fp_t MB = KB*1024.0;
const fp_t GB = MB*1024.0;

namespace qontra {

typedef std::vector<fp_t>   poly_t;

template <class T, class U, class V>
using TwoLevelMap = std::map<T, std::map<U, V>>;

// Two-level maps:
namespace tlm {

template <class T, class U, class V> void 
put(TwoLevelMap<T, U, V>& m, T x, U y, V z) {
    if (!m.count(x))    m[x] = std::map<U, V>();
    m[x][y] = z;
}

}   // tlm

//
// HELPER FUNCTIONS:
//
//  safe_create_directory creates folders for a directory (and any necessary parent folders).
//
//  is_subset_of checks if container A is a subset of container B.
//  xor_into add an element to a container if it does not exist, and removes it if it does.
//  to_vlw converts a stim::simd_bits or stim::simd_bits_range_ref object to a vlw_t

inline void
safe_create_directory(const std::filesystem::path& path) {
    if (path.string().size() == 0) return;
    if (!std::filesystem::exists(path)) {
        auto parent_path = path.parent_path();
        safe_create_directory(parent_path);
        std::filesystem::create_directory(path);
    }
}

// Checks if all elements in first arg are in second arg.
template <class T, class U> inline bool 
is_subset_of(T c1, U c2) {
    U tmp;
    std::set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(),
            std::inserter(tmp, tmp.begin()));
    return tmp.size() == c1.size();
}

template <class T> inline void
xor_into(std::set<T>& s, T x) {
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
}

template <typename T> int
number_of_common_elements(const std::vector<T>& arr1, const std::vector<T>& arr2) {
    std::set<T> ss(arr1.begin(), arr1.end());
    for (T x : arr2) {
        if (!ss.count(x)) return false;
        ss.erase(x);
    }
    return ss.size();
}

 // T must either be stim::simd_bits or stim::simd_bits_range_ref
template <class T> inline vlw_t
to_vlw(T x, uint64_t len) {
    uint n = (len >> 6)+1;
    vlw_t w(n);
    for (uint i = 0; i < n; i++)    w[i] = x.u64[i];
    return w;
}

inline int
vlw_compare(vlw_t a, vlw_t b) {
    int max_size = a.size() > b.size() ? a.size() : b.size();
    int min_size = a.size() < b.size() ? a.size() : b.size();
    bool a_is_larger = a.size() > b.size();
    for (int i = max_size-1; i >= 0; i--) {
        if (i >= min_size) {
            if (a_is_larger && a[i] > 0)        return 1;
            else if (!a_is_larger && b[i] > 0)  return -1;
        } else {
            if (a[i] > b[i])        return 1;
            else if (b[i] > a[i])   return -1;
        }
    }
    return 0;
}

}  // qontra

#endif
