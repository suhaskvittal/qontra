/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DEFS_h
#define DEFS_h

#include <filesystem>
#include <map>

#include <stdint.h>

typedef double fp_t;

#ifdef __APPLE__
typedef uint32_t uint;
#endif

typedef uint64_t addr_t;

const fp_t KB = 1024.0
const fp_t MB = KB*1024.0
const fp_t GB = MB*1024.0;

namespace qontra {

// Useful typedefs:
typedef std::vector<fp_t>   poly_t;

using TwoLevelMap<T, U, V> = std::map<T, std::map<U, V>>;

namespace tlm {

template <class T, class U, class V> void 
put(TwoLevelMap<T, U, V>& m, T x, U y, V z) {
    if (!m.count(x))    m[x] = std::map<U, V>();
    m[x][y] = z;
}

}   // tlm

void
safe_create_directory(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        auto parent_path = path.parent_path();
        safe_create_directory(parent_path);
        std::filesystem::create_directory(path);
    }
}

}  // qontra

#endif
