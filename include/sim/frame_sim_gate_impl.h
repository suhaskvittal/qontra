/*
 *  author: Suhas Vittal
 *  date:   14 December 2023
 * */

#ifndef FRAME_SIM_GATE_IMPL_h
#define FRAME_SIM_GATE_IMPL_h

// All functions below are just templated functions used
// by the FrameSimulator.

template <typename T> void
__h_gate(T& x, T& z, T& lock) {
    T tmp = z;
    z = (x & ~lock) | (z & lock);
    x = (tmp & ~lock) | (x & lock);
}

template <typename T> void
__x_gate(T& x, T& lock) {
    x ^= ~lock;
}

template <typename T> void
__z_gate(T& z, T& lock) {
    z ^= ~lock;
}

template <typename T, typename R> void
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
    x1 ^= l2 & rx1 & ~lock_or;
    z1 ^= ((z2 & ~l2) | (l2 & rz1)) & ~lock_or;
    x2 ^= ((x1 & ~l1) | (l1 & rx2)) & ~lock_or;
    z2 ^= l1 & rz2 & ~lock_or;
}

template <typename T> void
__liswap_gate(T& x1, T& x2, T& l2, T& lock) {
    // Swap in |02>/|11> subspace.
    auto s02_11 = ~x1 & l2;
    auto s11_02 = x1 & ~l2;

    x2 ^= s02_11 & ~lock;
    l2 ^= (s02_11 | s11_02) & ~lock;
    x1 ^= ~lock;
}

template <typename T, typename R> void
__measure(T& x, T& z, T& l, T& lock, R& rx, R& rz) {
    x = (((x & ~l) | (rx & l)) & ~lock) | (x & lock);
    z = (rz & ~lock) | (z & lock);
}

template <typename T, typename R> void
__reset(T& x, T& z, T& l, T& lock, R& rand) {
    x &= lock;
    l &= lock;
    z = (rand & ~lock) | (z & lock);
}


#endif  // FRAME_SIM_GATE_IMPL_h
