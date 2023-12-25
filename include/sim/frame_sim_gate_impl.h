/*
 *  author: Suhas Vittal
 *  date:   14 December 2023
 * */

#ifndef FRAME_SIM_GATE_IMPL_h
#define FRAME_SIM_GATE_IMPL_h

#include "stimext.h"

// All functions below are just templated functions used
// by the FrameSimulator.

template <class T> void
__h_gate(T& x, T& z, T& lock) {
    T tmp = z;
    z = andnot(lock, x) | (z & lock);
    x = andnot(lock, tmp) | (x & lock);
}

template <class T> void
__x_gate(T& x, T& lock) {
    x ^= ~lock;
}

template <class T> void
__z_gate(T& z, T& lock) {
    z ^= ~lock;
}

template <class T, class R> void
__cx_gate(T& x1,
            T& z1,
            T& l1,
            T& x2,
            T& z2,
            T& l2,
            T& lock_or,
            R& rx1,
            R& rz1,
            R& rx2,
            R& rz2)
{
    x1 ^= andnot(lock_or, l2 & rx1);
    z1 ^= andnot(lock_or, ((z2 & ~l2) | (l2 & rz1)));
    x2 ^= andnot(lock_or, ((x1 & ~l1) | (l1 & rx2)));
    z2 ^= andnot(lock_or, l1 & rz2);
}

template <class T> void
__liswap_gate(T& x1, T& x2, T& l2, T& lock) {
    // Swap in |02>/|11> subspace.
    auto s02_11 = ~x1 & l2;
    auto s11_02 = x1 & ~l2;

    x2 ^= andnot(lock, s02_11);
    l2 ^= andnot(lock, (s02_11 | s11_02));
    x1 ^= ~lock;
}

template <class T, class R> void
__measure(T& x, T& z, T& l, T& lock, R& rx, R& rz) {
    x = andnot(lock, (x & ~l) | (rx & l)) | (x & lock);
    z = andnot(lock, rz) | (z & lock);
}

template <class T, class R> void
__reset(T& x, T& z, T& l, T& lock, R& rand) {
    x &= lock;
    l &= lock;
    z = andnot(lock, rand) | (z & lock);
}


#endif  // FRAME_SIM_GATE_IMPL_h
