/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#ifndef VTILS_FIXED_STRING_h
#define VTILS_FIXED_STRING_h

// Useful for string templates.

#include <algorithm>

namespace vtils {

template <size_t N>
struct fixed_string {
    constexpr fixed_string(const char (&s)[N]) {
        std::copy_n(s, N, value);
    }

    char value[N];
};

}   // vtils

#endif  // VTILS_FIXED_STRING_h
