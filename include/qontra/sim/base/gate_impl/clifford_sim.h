/*
 *  author: Suhas Vittal
 *  date:   29 January 2024
 * */

#ifndef CLIFFORD_SIM_GATE_IMPL_h
#define CLIFFORD_SIM_GATE_IMPL_h

namespace qontra {

template <class T> inline void
__h_gate(T& r, T& x, T& z, T& lock) {
    T tmp = z;
    r ^= andnot(lock, x & z);
    z = (lock & z) | andnot(lock, x);
    x = (lock & x) | andnot(lock, tmp);
}

template <class T> inline void
__x_gate(T& r, T& z, T& lock) {
    r ^= andnot(lock, z);
}

template <class T> inline void
__z_gate(T& r, T& x, T& lock) {
    r ^= andnot(lock, x);
}

template <class T> inline void
__s_gate(T& r, T& x, T& z, T& lock) {
    r ^= andnot(lock, x & z);
    z ^= andnot(lock, x);
}

template <class T> inline void
__cx_gate(T& r, T& x1, T& z1, T& x2, T& z2, T& lock) {
    r ^= andnot(lock, andnot(x2 ^ z1, (x1 & z2)));
    x2 ^= andnot(lock, x1);
    z1 ^= andnot(lock, z2);
}

template <class T, class U> inline void
__rowsum_arith(T& x1, T& z1, T& x2, T& z2, T& s1, T& s2, U p) {
    // Truth Table:
    // |----|----|-------|
    // | 1  | 2  |   g   |
    // |----|----|-------|
    // | 00 | xx |   0   |
    // | 01 | 0x |   0   |
    // | 01 | 10 |   1   |  --> mag
    // | 01 | 11 |  -1   |  --> sgn
    // | 10 | x0 |   0   |
    // | 10 | 01 |  -1   |  --> sgn
    // | 10 | 11 |   1   |  --> mag
    // | 11 | 00 |   0   |
    // | 11 | 01 |   1   |  --> mag
    // | 11 | 10 |  -1   |  --> sgn
    // | 11 | 11 |   0   |
    // |----|----|-------|
    //
    // I'll be honest, I completely forgot what the hell I did here.
    // But whatever me wrote this was a genius because it works!
    //  --> Well, it'd be awkward if it does not actually work...
    //
    // Basic idea: mag and sgn track the value of the rowsum, and we can 
    // implement Z_4 arithmetic using bitwise operations.
    auto mag = (andnot(x1, z1) & x2)
                | (andnot(z1, x1) & z2)
                | (x1 & z1 & (x2 ^ z2));
    auto sgn = (andnot(x1, z1) & x2 & z2)
                | (andnot(z1, x1) & x2 & z2)
                | (x1 & z1 & andnot(z2, x2));
    s2 ^= p & ((s1 & mag) ^ sgn);
    s1 ^= p & mag;
    x2 ^= p & x1;
    z2 ^= p & z1;
}

}   // qontra

#endif  // CLIFFORD_SIM_GATE_IMPL_h
