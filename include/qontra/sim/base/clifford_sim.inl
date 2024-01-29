/*
 *  author: Suhas Vittal
 *  date:   29 January 2024
 * */

namespace qontra {

inline void
CliffordSimulator::reset_sim(void) {
    StateSimulator::reset_sim();

    x_table.clear();
    z_table.clear();
    r_table.clear();
    init_tables();
}

inline void
CliffordSimulator::eDP1(uint64_t q, uint64_t t) {
    auto p = rng() & 3;
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        r_table[i][t] ^= (z_table[k][t] & (p & 1)) ^ (x_table[k][t] & (p & 2));
    }
}

inline void
CliffordSimulator::eX(uint64_t q, uint64_t t) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        r_table[i][t] ^= z_table[k][t];
    }
}

inline void
CliffordSimulator::eY(uint64_t q, uint64_t t) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        r_table[i][t] ^= z_table[k][t] ^ x_table[k][t];
    }
}

inline void
CliffordSimulator::eZ(uint64_t q, uint64_t t) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        r_table[i][t] ^= x_table[k][t];
    }
}

inline void
CliffordSimulator::eL(uint64_t q, uint64_t t) { }

inline void
CliffordSimulator::eDP2(uint64_t q1, uint64_t q2, uint64_t t) {
    auto p = rng() & 15;
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        uint64_t k1 = get_index(i, q1),
                 k2 = get_index(i, q2);
        r_table[i][t] ^= (z_table[k1][t] & (p & 1)) ^ (x_table[k1][t] & (p & 2));
        r_table[i][t] ^= (z_table[k2][t] & (p & 4)) ^ (x_table[k2][t] & (p & 8));
    }
}

inline void
CliffordSimulator::eLI(uint64_t q1, uint64_t q2, uint64_t t) { }

inline void
CliffordSimulator::eLT(uint64_t q1, uint64_t q2, uint64_t t) { }

inline void
CliffordSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = x_table;
    z_table_cpy = z_table;
    r_table_cpy = r_table;
}

inline void
CliffordSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    StateSimulator::rollback_where(pred);
    copy_where(x_table_cpy, x_table, pred);
    copy_where(z_table_cpy, z_table, pred);
    copy_where(r_table_cpy, r_table, pred);
}

inline uint64_t
CliffordSimulator::get_index(uint64_t i, uint64_t q) {
    return i*n_qubits + q;
}

inline void
CliffordSimulator::clear_scratch_space() {
    // Clear scratch space
    for (uint64_t i = 2*n_qubits*n_qubits; i < x_width; i++)    x_table[i].clear();
    for (uint64_t i = 2*n_qubits*n_qubits; i < z_width; i++)    z_table[i].clear();
    for (uint64_t i = 2*n_qubits; i < r_width; i++)             r_table[i].clear();
}

inline void
CliffordSimulator::clear_row(uint64_t r) {
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t k = get_index(r, j);
        x_table[k].clear();
        z_table[k].clear();
    }
    r_table[r].clear();
}

inline void
CliffordSimulator::clear_row_for_trial(uint64_t r, int64_t tr) {
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t k = get_index(r, j);
        x_table[k][tr] = 0;
        z_table[k][tr] = 0;
    }
    r_table[r][tr] = 0;
}

inline void
CliffordSimulator::clear_row_where(uint64_t r, stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    pred.invert_bits();
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t k = get_index(r, j);
        x_table[k] &= pred;
        z_table[k] &= pred;
    }
    r_table[r] &= pred;
    pred.invert_bits();
}

inline void
CliffordSimulator::swap_rows_for_trial(uint64_t r1, uint64_t r2, int64_t tr) {
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t k1 = get_index(r1, j),
                 k2 = get_index(r2, j);
        x_table[k1][tr].swap_with(x_table[k2][tr]);
        z_table[k1][tr].swap_with(z_table[k2][tr]);
    }
    r_table[r1][tr].swap_with(r_table[r2][tr]);
}

inline void
CliffordSimulator::swap_rows_where(uint64_t r1, uint64_t r2, stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t k1 = get_index(r1, j),
                 k2 = get_index(r2, j);
        x_table[k1].for_each_word(
                x_table[k2], z_table[k1], z_table[k2], pred,
            [] (auto& x1, auto& x2, auto& z1, auto& z2, auto& p)
            {
                auto tmpx = x1,
                     tmpz = z1;
                x1 = (p & x2) | andnot(p, x1);
                x2 = (p & tmpx) | andnot(p, x2);
                z1 = (p & z2) | andnot(p, z1);
                z2 = (p & tmpz) | andnot(p, z2);
            });
    }
    r_table[r1].for_each_word(r_table[r2], pred,
            [] (auto& r1, auto& r2, auto& p)
            {
                auto tmp = r1;
                r1 = (p & r2) | andnot(p, r1);
                r2 = (p & tmp) | andnot(p, r2);
            });
}

inline void
CliffordSimulator::init_tables() {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (uint64_t j = 0; j < n_qubits; j++) {
            uint64_t k = get_index(i, j);
            if (i < n_qubits && i == j) {
                x_table[k].invert_bits();
            } else if ((i-n_qubits) == j) {
                z_table[k].invert_bits();
            }
        }
    }
}

}   // qontra
