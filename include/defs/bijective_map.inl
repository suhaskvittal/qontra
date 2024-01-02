/*
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#include <type_traits>

#include <stdlib.h>

template <class T, class U> inline BijectiveMap<T, U>&
BijectiveMap<T, U>::operator=(const BijectiveMap<T, U>& other) {
    forward_map = other.forward_map;
    reverse_map = other.reverse_map;
    return *this;
}

template <class T, class U> inline BijectiveMap<T, U>
BijectiveMap<T, U>::reverse() {
    return from_maps(reverse_map, forward_map);
}

template <class T, class U> inline bool
BijectiveMap<T, U>::f_put(T key, U value) {
    if (f_count(key) && r_count(value) && f_at(key) != value) {
        // Then, we have a situation:
        //  f(key) = X, r(X) = key
        //  r(value) = Y, f(Y) = value
        // Post-update:
        //  f(key) = value, r(X) = key
        //  r(value) = key, f(Y) = value
        if (assert_bijectivity) exit(1);
        return false;
    } else if (f_count(key) && !r_count(value)) {
        reverse_map.erase(f_at(key));   // Delete existing entry in reverse map.
    } else if (!f_count(key) && r_count(value)) {
        forward_map.erase(r_at(value));
    }
    forward_map[key] = value;
    reverse_map[value] = key;
    return true;
}

template <class T, class U> inline bool
BijectiveMap<T, U>::r_put(U key, T value) {
    return f_put(value, key);
}

template <class T, class U> inline void
BijectiveMap<T, U>::f_erase(T key) {
    U value = f_at(key);
    forward_map.erase(key);
    reverse_map.erase(value);
}

template <class T, class U> inline void
BijectiveMap<T, U>::r_erase(U key) {
    T value = r_at(key);
    reverse_map.erase(key);
    forward_map.erase(value);
}

template <class T, class U> inline void
BijectiveMap<T, U>::f_swap(T k1, T k2) {
    U v1 = f_at(k1),
      v2 = f_at(k2);
    std::swap(forward_map[k1], forward_map[k2]);
    std::swap(reverse_map[v1], reverse_map[v2]);
}

template <class T, class U> inline void
BijectiveMap<T, U>::r_swap(U k1, U k2) {
    T v1 = r_at(k1),
      v2 = r_at(k2);
    std::swap(reverse_map[k1], reverse_map[k2]);
    std::swap(forward_map[v1], forward_map[v2]);
}
