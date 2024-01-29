/*
 *  author: Suhas Vittal
 * */

namespace qontra {

inline void
FrameSimulator::reset_sim() {
    StateSimulator::reset_sim();

    x_table.clear();
    z_table.clear();
    leak_table.clear();
}

inline void
FrameSimulator::eDP1(uint64_t q, uint64_t t) {
    int p = rng() & 3;
    x_table[q][t] ^= (bool)(p & 1);
    z_table[q][t] ^= (bool)(p & 2);
}

inline void
FrameSimulator::eX(uint64_t q, uint64_t t) {
    x_table[q][t] ^= 1;
}

inline void
FrameSimulator::eY(uint64_t q, uint64_t t) {
    x_table[q][t] ^= 1;
    z_table[q][t] ^= 1;
}

inline void
FrameSimulator::eZ(uint64_t q, uint64_t t) {
    z_table[q][t] ^= 1;
}

inline void
FrameSimulator::eL(uint64_t q, uint64_t t) {
    leak_table[q][t] ^= 1;
}

inline void
FrameSimulator::eDP2(uint64_t q1, uint64_t q2, uint64_t t) {
    int p = rng() & 15;
    x_table[q1][t] ^= (bool)(p & 1);
    z_table[q1][t] ^= (bool)(p & 2);
    x_table[q2][t] ^= (bool)(p & 4);
    z_table[q2][t] ^= (bool)(p & 8);
}

inline void
FrameSimulator::eLI(uint64_t q1, uint64_t q2, uint64_t t) {
    int p = rng() % 3;
    bool c1 = (p == 0) || (p == 2);
    bool c2 = (p == 1) || (p == 2);
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
}

#define QONTRA_USE_MESSY_LEAKAGE_TRANSPORT

inline void 
FrameSimulator::eLT(uint64_t q1, uint64_t q2, uint64_t t) {
#ifdef QONTRA_USE_MESSY_LEAKAGE_TRANSPORT
    bool c1 = leak_table[q2][t] & !leak_table[q1][t];
    bool c2 = leak_table[q1][t] & !leak_table[q2][t];
#else
    bool c1 = leak_table[q1][t] ^ leak_table[q2][t];
    bool c2 = c1;
#endif
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
}

inline void
FrameSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(x_table);
    z_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(z_table);
    leak_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(leak_table);
}

inline void
FrameSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    StateSimulator::rollback_where(pred);
    for (size_t i = 0; i < n_qubits; i++) {
        copy_where(x_table_cpy[i], x_table[i], pred);
        copy_where(z_table_cpy[i], z_table[i], pred);
        copy_where(leak_table_cpy[i], leak_table[i], pred);
    }
}

}   // qontra
