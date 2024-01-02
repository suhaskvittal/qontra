/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef DEFS_SET_ALGEBRA_h
#define DEFS_SET_ALGEBRA_h

#include <algorithm>
#include <set>

// Useful operator overloads:
namespace std {

template <class T> inline set<T>
operator+(set<T> s1, set<T> s2) {
    std::set<T> s3;
    std::set_union(s1.begin(),
                    s1.end(),
                    s2.begin(),
                    s2.end(),
                    std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline set<T>
operator-(set<T> s1, set<T> s2) {
    std::set<T> s3;
    std::set_difference(s1.begin(),
                        s1.end(),
                        s2.begin(),
                        s2.end(),
                        std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline set<T>
operator*(set<T> s1, set<T> s2) {
    std::set<T> s3;
    std::set_intersection(s1.begin(),
                            s1.end(),
                            s2.begin(),
                            s2.end(),
                            std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline set<T>
operator^(set<T> s1, set<T> s2) {
    std::set<T> s3;
    std::set_symmetric_difference(s1.begin(),
                                    s1.end(),
                                    s2.begin(),
                                    s2.end(),
                                    std::inserter(s3, s3.begin()));
    return s3;
}

template <class T> inline set<T>&
operator+=(set<T>& s1, set<T> s2) { s1 = s1 + s2; return s1; }

template <class T> inline set<T>&
operator-=(set<T>& s1, set<T> s2) { s1 = s1 - s2; return s1; }

template <class T> inline set<T>&
operator*=(set<T>& s1, set<T> s2) { s1 = s1 * s2; return s1; }

template <class T> inline set<T>&
operator^=(set<T>& s1, set<T> s2) { s1 = s1 ^ s2; return s1; }

template <class T> inline set<T>&
operator+=(set<T>& s, T x) { s.insert(x); return s; }

template <class T> inline set<T>&
operator-=(set<T>& s, T x) { s.erase(x); return s; }

template <class T> inline set<T>&
operator^=(set<T>& s, T x) { 
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
    return s;
}

}   // std

#endif  // DEFS_SET_ALGEBRA_h
