/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef DEFS_SET_ALGEBRA_h
#define DEFS_SET_ALGEBRA_h

#include <algorithm>
#include <set>

template <class T> inline std::set<T>
set_union(std::set<T> s1, std::set<T> s2) {
    std::set<T> s3;
    std::set_union(s1.begin(),
                    s1.end(),
                    s2.begin(),
                    s2.end(),
                    std::inserter(s3, s3.begin()));
}

template <class T> inline std::set<T>
set_intersect(std::set<T> s1, std::set<T> s2) {
    std::set<T> s3;
    std::set_intersection(s1.begin(),
                            s1.end(),
                            s2.begin(),
                            s2.end(),
                            std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline std::set<T>
set_difference(std::set<T> s1, std::set<T> s2) {
    std::set<T> s3;
    std::set_difference(s1.begin(),
                        s1.end(),
                        s2.begin(),
                        s2.end(),
                        std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline std::set<T>
set_symmetric_difference(std::set<T> s1, std::set<T> s2) {
    std::set<T> s3;
    std::set_symmetric_difference(s1.begin(),
                                    s1.end(),
                                    s2.begin(),
                                    s2.end(),
                                    std::inserter(s3, s3.begin()));
    return s3;
}

// Useful operator overloads:
namespace std {

template <class T> inline std::set<T>
operator+(std::set<T> s1 std::set<T> s2) { return set_union(s1, s2); }

template <class T> inline std::set<T>
operator-(std::set<T> s1 std::set<T> s2) { return set_difference(s1, s2); }

template <class T> inline std::set<T>
operator*(std::set<T> s1 std::set<T> s2) { return set_intersection(s1, s2); }

template <class T> inline std::set<T>
operator^(std::set<T> s1 std::set<T> s2) { return set_symmetric_difference(s1, s2); }

template <class T> inline std::set<T>&
operator+=(std::set<T>& s1, std::set<T> s2) { s1 = s1 + s2; return s1; }

template <class T> inline std::set<T>&
operator-=(std::set<T>& s1, std::set<T> s2) { s1 = s1 - s2; return s1; }

template <class T> inline std::set<T>&
operator*=(std::set<T>& s1, std::set<T> s2) { s1 = s1 * s2; return s1; }

template <class T> inline std::set<T>&
operator^=(std::set<T>& s1, std::set<T> s2) { s1 = s1 ^ s2; return s1; }

template <class T> inline std::set<T>&
operator+=(std::set<T>& s, T x) { s.insert(x); return s; }

template <class T> inline std::set<T>&
operator-=(std::set<T>& s, T x) { s.erase(x); return s; }

template <class T> inline std::set<T>&
operator^=(std::set<T>& s, T x) { 
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
    return s;
}

}   // std

#endif  // DEFS_SET_ALGEBRA_h
