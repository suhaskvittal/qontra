/*
 *  author: Suhas Vittal
 *  date:   18 June 2023
 * */

#include "sim/clifford_sim.h"

namespace qontra {

#define L_I     StateSimulator::label_1q_t::I
#define L_X     StateSimulator::label_1q_t::X
#define L_Y     StateSimulator::label_1q_t::Y
#define L_Z     StateSimulator::label_1q_t::Z
#define L_L     StateSimulator::label_1q_t::L

static const StateSimulator::label_1q_t IXYZ[4] = { L_I, L_X, L_Y, L_Z };
static const StateSimulator::label_1q_t IX[2] = { L_I, L_X };
static const StateSimulator::label_1q_t IY[2] = { L_I, L_Y };
static const StateSimulator::label_1q_t IZ[2] = { L_I, L_Z };
static const StateSimulator::label_1q_t IL[2] = { L_I, L_L };

void
CliffordSimulator::H(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;    // Note that x_table and z_table are
                                        // flattened 2D arrays.
            r_table[i].for_each_word(
                x_table[k], z_table[k], leak_table[j], lock_table[j],
            [&](auto& r, auto& x, auto& z, auto& l, auto& lock)
            {
                r ^= x & z & ~l & ~lock;
                stim::simd_word tmp = z;
                z = (z & (l | lock)) | (x & ~(l | lock));
                x = (x & (l | lock)) | (tmp & ~(l | lock));
            });
        }
    }
}

void
CliffordSimulator::X(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            // Only flip where not leaked.
            r_table[i].for_each_word(
                z_table[k], leak_table[j], lock_table[j],
            [&](auto& r, auto& z, auto& l, auto& lock)
            {
                r ^= z & ~l & ~lock;
            });
        }
    }
}

void
CliffordSimulator::Z(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            r_table[i].for_each_word(
                x_table[k], leak_table[j], lock_table[j],
            [&](auto& r, auto& x, auto& l, auto& lock)
            {
                r ^= x & ~l & ~lock;
            });
        }
    }
}

void
CliffordSimulator::S(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            r_table[i].for_each_word(
                    x_table[k], z_table[k], leak_table[j], lock_table[j],
            [&](auto& r, auto& x, auto& z, auto& l, auto& lock)
            {
                r ^= x & z & ~l & ~lock;
                z ^= x & ~l & ~lock;
            });
        }
    }
}

void
CliffordSimulator::CX(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint i = 0; i < operands.size(); i += 2) {
            uint j1 = operands[i];
            uint j2 = operands[i+1];
            uint k1 = i*n_qubits + j1;
            uint k2 = i*n_qubits + j2;

            stim::simd_bits lock_both(lock_table[j1]);
            lock_both |= lock_table[j2];

            r_table[i].for_each_word(
                x_table[k1], z_table[k1],
                x_table[k2], z_table[k2],
                leak_table[j1], leak_table[j2],
                lock_both,
            [&](auto& r, auto& x1, auto& z1, 
                auto& x2, auto& z2, 
                auto& l1, auto& l2,
                auto& lock)
            {
                stim::simd_word rx1(rng(), rng());
                stim::simd_word rz1(rng(), rng());
                stim::simd_word rx2(rng(), rng());
                stim::simd_word rz2(rng(), rng());

                // Only update phase if the qubits are not leaked.
                r ^= ((x1 & z2) & ~(x2 ^ z1) & ~(l1 | l2)) & ~lock;
                x2 ^= ~(l1 | l2 | lock) & x1;
                z1 ^= ~(l1 | l2 | lock) & z2;
                // Apply depolarizing error if they are leaked.
                r ^= ((z1 & rx1 & l2)
                        & (x1 & rz1 & l2)
                        & (z2 & rx2 & l1)
                        & (x2 & rz2 & l2)) & ~lock;
            });
        }
    }
}

void
CliffordSimulator::M(
        std::vector<uint> operands,
        fp_t m1w0,
        fp_t m0w1,
        int record)
{
    for (uint j : operands) {
        // Clear scratch space
        for (uint i = 2*n_qubits*n_qubits; i < x_width; i++)    x_table[i].clear();
        for (uint i = 2*n_qubits*n_qubits; i < z_width; i++)    z_table[i].clear();
        for (uint i = 2*n_qubits; i < r_width; i++)             r_table[i].clear();

        // Figure out if xij is nonzero anywhere (just need to OR all rows).
        for (uint i = n_qubits; i < 2*n_qubits; i++) {
            uint k = n_qubits*i + j;
            x_table[x_width-1] |= x_table[k];
        }
        // If the qubit was leaked, also set it to 0: we can quickly determine
        // the random result here.
        leak_table[j].invert_bits();
        x_table[x_width-1] &= leak_table[j];
        leak_table[j].invert_bits();

        // If the qubit is locked in a trial, also set it to 0. We will just
        // record the measurement for the locked trial as 0.
        lock_table[j].invert_bits();
        x_table[x_width-1] &= lock_table[j];
        
        // Perform determinate measurement. We want to condition this
        // to happen wherever x_table[x_width-1] = 0
        // 
        // Clear out row 2*n_qubits+1
        for (uint jj = 0; jj < n_qubits; jj++) {
            x_table[2*n_qubits*n_qubits + jj].clear();  // Remember X and Z tables are
            z_table[2*n_qubits*n_qubits + jj].clear();  // flattened 2d arrays.
        }
        r_table[2*n_qubits].clear();
        for (uint i = 0; i < n_qubits; i++) {
            // Predicate the rowsum where xij = 1 and j is unlocked.
            uint k = n_qubits*i + j;
            stim::simd_bits pred(x_table[k]);
            pred &= lock_table[j];
            browsum(2*n_qubits, i+n_qubits, true, pred);
        }
        lock_table[j].invert_bits();
        // Measurement outcome is r_table[2*n_qubits]
        // Condition results based on x_table[x_width-1] = 0
        // Also, if any of these qubits were leaked, set the measurement
        // output to a random value.
        r_table[2*n_qubits].for_each_word(
                x_table[x_width-1], leak_table[j], lock_table[j],
        [&] (auto& r, auto& x, auto& l, auto& lock)
        {
            stim::simd_word rand(rng(), rng());
            r = ((r & ~l) | (rand & l)) & ~(x | lock);
        });
        // Now, perform any indeterminate measurements.
        // This takes more time as we must handle each case separately.
        for (uint64_t t = 0; t < shots; t++) {
            if (!x_table[x_width-1][t]) continue;
            // Find first i where xia = 1 in i = n+1 to 2n. Call this ii.
            uint i;
            for (uint i = n_qubits; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (x_table[k][t]) {
                    i = i;
                    break;
                }
            }
            // First perform rowsums.
            for (uint i = 0; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (i != i && x_table[k][t]) {
                    rowsum(i, i, t);
                }
            }
            // Second, swap (ii-n)-th row with ii-th row and clear
            // ii-th row. Set rii to random value.
            for (uint jj = 0; jj < n_qubits; jj++) {
                x_table[(i-n_qubits)*n_qubits+jj][t] = x_table[i*n_qubits+jj][t];
                z_table[(i-n_qubits)*n_qubits+jj][t] = z_table[i*n_qubits+jj][t];
                x_table[i*n_qubits+jj][t] = 0;
                z_table[i*n_qubits+jj][t] = 0;
            }
            r_table[i-n_qubits][t] = r_table[i][t];
            r_table[i][t] = rng() & 1;
            z_table[i*n_qubits+j][t] = 1;
            // For simplicity, copy the rii value (the measurement outcome)
            // to the scratch space (row 2n+1) in r.
            r_table[2*n_qubits][t] = r_table[i][t];
        }

        if (record >= 0) {
            record_table[record].for_each_word(r_table[2*n_qubits], lock_table[j],
                    [&] (auto& rec, auto& r, auto& lock)
                    {
                        rec |= r & ~lock;
                    });
            r_table[2*n_qubits].clear();
            error_channel_m(record, m1w0, m0w1, lock_table[j]);
            record++;
        }
    }
}

void
CliffordSimulator::R(std::vector<uint> operands) {
    // Remove leakage.
    for (uint j : operands) {
        leak_table[j].clear();
    }
    // Implement as a measure + X gate.
    stim::simd_bits record_0_cpy(record_table[0]);
    for (uint j : operands) {
        M({j}, 0);
        for (uint i = 0; i < 2*n_qubits; i++) {
            uint k = i*n_qubits + j;
            // Rec should be 0 whenever j is locked, so 
            // the outcome should not change.
            r_table[i].for_each_word(
                z_table[k], record_table[0],
            [&](auto& r, auto& z, auto& rec) {
                r ^= z & rec;
            });
        }
    }
    record_table[0].swap_with(record_0_cpy);
}

StateSimulator::label_1q_t
CliffordSimulator::eDP1(uint q, uint64_t t) {
    auto p = rng() & 3;
    if (leak_table[q][t])   p = 0;
    for (uint i = 0; i < 2*n_qubits; i++) {
        uint k = n_qubits*i + q;
        r_table[i][t] ^= (z_table[k][t] & (p & 1)) ^ (x_table[k][t] & (p & 2));
    }
    return IXYZ[p];
}

StateSimulator::label_1q_t
CliffordSimulator::eX(uint q, uint64_t t) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        uint k = n_qubits*i + q;
        r_table[i][t] ^= z_table[k][t] & ~leak_table[q][t];
    }
    return IX[1 - leak_table[q][t]];
}

StateSimulator::label_1q_t
CliffordSimulator::eY(uint q, uint64_t t) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        uint k = n_qubits*i + q;
        r_table[i][t] ^= (x_table[k][t] ^ z_table[k][t]) & ~leak_table[q][t];
    }
    return IY[1 - leak_table[q][t]];
}

StateSimulator::label_1q_t
CliffordSimulator::eZ(uint q, uint64_t t) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        uint k = n_qubits*i + q;
        r_table[i][t] ^= x_table[k][t] & ~leak_table[q][t];
    }
    return IZ[1 - leak_table[q][t]];
}

StateSimulator::label_1q_t
CliffordSimulator::eL(uint q, uint64_t t) {
    leak_table[q][t] ^= 1;
    return L_L;
}

StateSimulator::label_2q_t
CliffordSimulator::eDP2(uint q1, uint q2, uint64_t t) {
    auto p = rng() & 15;
    if (leak_table[q1][t])  p &= 0xc;   // Clear out bottom 2 bits
    if (leak_table[q2][t])  p &= 0x3;   // Clear out upper two bits.
    for (uint i = 0; i < 2*n_qubits; i++) {
        uint k1 = n_qubits*i + q1;
        uint k2 = n_qubits*i + q2;
        r_table[i][t] ^= (z_table[k1][t] & (p & 1)) ^ (x_table[k1][t] & (p & 2));
        r_table[i][t] ^= (z_table[k2][t] & (p & 4)) ^ (x_table[k2][t] & (p & 8));
    }
    return std::make_pair(IXYZ[p & 3], IXYZ[(p >> 2) & 3]);
}

StateSimulator::label_2q_t 
CliffordSimulator::eLI(uint q1, uint q2, uint64_t t) {
    auto p = rng() % 3;
    bool c1 = (p == 0) || (p == 2);
    bool c2 = (p == 1) || (p == 2);
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
    return std::make_pair(IL[c1], IL[c2]);
}

StateSimulator::label_2q_t 
CliffordSimulator::eLT(uint q1, uint q2, uint64_t t) {
    bool c1 = leak_table[q2][t] & ~leak_table[q1][t];
    bool c2 = leak_table[q1][t] & ~leak_table[q2][t];
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
    return std::make_pair(IL[c1], IL[c2]);
}

void
CliffordSimulator::init_tables() {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j = 0; j < n_qubits; j++) {
            uint k = i*n_qubits + j;
            if (i < n_qubits && i == j) {
                x_table[k].invert_bits();
            } else if ((i-n_qubits) == j) {
                z_table[k].invert_bits();
            }
        }
    }
}

void
CliffordSimulator::rowsum(uint h, uint i, uint64_t t) {
    auto rs1 = 0;
    auto rs2 = r_table[h][t] ^ r_table[i][t];
    for (uint j = 0; j < n_qubits; j++) {
        uint kh = n_qubits*h + j;
        uint ki = n_qubits*i + j;

        auto x1 = x_table[ki][t];
        auto x2 = x_table[kh][t];
        auto z1 = z_table[ki][t];
        auto z2 = z_table[kh][t];
        
        auto g1 = (~x1 & z1 & x2)
                    | (x1 & ~z1 & z2) | (x1 & z1 & (x2 ^ z2));
        auto g2 = (~x1 & z1 & x2 & x2) 
                    | (x1 & ~z1 & x2 & z2)
                    | (x1 & z1 & x2 & ~z2);
        rs2 ^= (rs1 & g1) ^ g2;
        rs1 ^= g1;
        // Update X and Z
        x_table[kh][t] ^= x1;
        z_table[kh][t] ^= z1;
    }
    r_table[h][t] = rs2;
}

void
CliffordSimulator::browsum(uint h, uint i, bool use_pred, 
        stim::simd_bits_range_ref pred) 
{
    // 2n + 1 and 2n + 2 correspond to the 1s and 2s place in a 
    // 2-bit number.
    //
    // To perform the 2rh + 2ri computation, we can just
    // XOR the bits going to the 2s place (as this is mod 4).
    stim::simd_word en_p(0, 0);
    if (use_pred)   en_p = ~en_p;

    r_table[2*n_qubits+2] ^= r_table[h];
    r_table[2*n_qubits+2] ^= r_table[i];
    for (uint j = 0; j < n_qubits; j++) {
        uint kh = n_qubits * h + j;
        uint ki = n_qubits * i + j;
        x_table[ki].for_each_word(
            z_table[ki], x_table[kh], z_table[kh],
            r_table[2*n_qubits+1], r_table[2*n_qubits+2],
            pred,
        [&](auto& x1, auto& z1, auto& x2, auto& z2, auto& scr1, auto& scr2, auto& p)
        {
            // Table:
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
            // I'm assuming that gsum_mag and gsum_sgn track the value of
            // the rowsum, and we can implement Z_4 arithmetic using bitwise
            // operations. I'm sure it's more obvious once you sit down with
            // pen and paper.

            stim::simd_word gmag(0, 0);
            stim::simd_word gsgn(0, 0);
            
            gmag ^= (~x1 & z1 & x2)
                        | (x1 & ~z1 & z2) | (x1 & z1 & (x2 ^ z2));
            gsgn ^= (~x1 & z1 & x2 & x2) 
                        | (x1 & ~z1 & x2 & z2)
                        | (x1 & z1 & x2 & ~z2);
            // Perform the modular addition and save the results in
            // the scratch storage
            scr2 ^= (en_p & p) & ((scr1 & gmag) ^ gsgn);
            scr1 ^= (en_p & p) & gmag;
            // Update x2 and z2
            x2 ^= (en_p & p) & x1;
            z2 ^= (en_p & p) & z1;
        });
    }   // After all the iterations, row 2n+1 and 2n+2 should have 
        // the value of the rowsum.

    // Now, in locations where row 2n+2 == 0, we need to set rh = 0. Otherwise,
    // we set it to 1. This is as simple as just moving row 2n+2 to rh.
    r_table[h].for_each_word(r_table[2*n_qubits+2], pred,
    [&] (auto& rh, auto& scr, auto& p) {
        rh = ((en_p & p) & scr) | (~(en_p & p) & rh);
    });
    // Clear the scratch storage used.
    r_table[2*n_qubits+1].clear();
    r_table[2*n_qubits+2].clear();
}

void
CliffordSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = stim::simd_bit_table(x_table);
    z_table_cpy = stim::simd_bit_table(z_table);
    r_table_cpy = stim::simd_bit_table(r_table);
    leak_table_cpy = stim::simd_bit_table(leak_table);
}

void
CliffordSimulator::rollback_where(stim::simd_bits_range_ref pred) {
    StateSimulator::rollback_where(pred);
    for (uint i = 0; i < x_width; i++)  copy_where(x_table_cpy[i], x_table[i], pred);
    for (uint i = 0; i < z_width; i++)  copy_where(z_table_cpy[i], z_table[i], pred);
    for (uint i = 0; i < r_width; i++)  copy_where(r_table_cpy[i], r_table[i], pred);
    for (uint i = 0; i < leak_width; i++) {
        copy_where(leak_table_cpy[i], leak_table[i], pred);
    }
}

}   // qontra
