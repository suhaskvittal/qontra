/* author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DEFS_h
#define DEFS_h

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

typedef uint64_t addr_t;

const fp_t KB = 1024.0;
const fp_t MB = KB*1024.0;
const fp_t GB = MB*1024.0;

namespace qontra {

// Useful typedefs:
typedef std::vector<fp_t>   poly_t;

template <class T, class U, class V>
using TwoLevelMap = std::map<T, std::map<U, V>>;

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
template <class T, class U>
inline bool is_subset_of(T c1, U c2) {
    U tmp;
    std::set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(),
            std::inserter(tmp, tmp.begin()));
    return tmp.size() == c1.size();
}

}  // qontra

#endif
