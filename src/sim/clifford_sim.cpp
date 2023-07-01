/*
 *  author: Suhas Vittal
 *  date:   18 June 2023
 * */

#include "sim/clifford_sim.h"

namespace qontra {

void
CliffordSimulator::H(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;    // Note that x_table and z_table are
                                        // flattened 2D arrays.
            r_table[i].for_each_word(
                x_table[k], z_table[k], leak_table[j],
            [&](auto& r, auto& x, auto& z, auto& l)
            {
                r ^= x & z & ~l;
                stim::simd_word tmp = z;
                z = (z & l) | (x & ~l);
                x = (x & l) | (tmp & ~l);
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
                z_table[k], leak_table[j],
            [&](auto& r, auto& z, auto& l)
            {
                r ^= z & ~l;
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
                x_table[k], leak_table[j],
            [&](auto& r, auto& x, auto& l)
            {
                r ^= x & ~l;
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
                    x_table[k], z_table[k], leak_table[j],
            [&](auto& r, auto& x, auto& z, auto& l)
            {
                r ^= x & z & ~l;
                z ^= x & ~l;
            });
        }
    }
}

void
CliffordSimulator::CX(std::vector<uint> operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint ii = 0; ii < operands.size(); ii += 2) {
            uint j1 = operands[ii];
            uint j2 = operands[ii+1];
            uint k1 = i*n_qubits + j1;
            uint k2 = i*n_qubits + j2;

            r_table[i].for_each_word(
                x_table[k1], z_table[k1],
                x_table[k2], z_table[k2],
                leak_table[j1], leak_table[j2],
            [&](auto& r, auto& x1, auto& z1, auto& x2, auto& z2, auto& l1, auto& l2)
            {
                stim::simd_word rx1(rng(), rng());
                stim::simd_word rz1(rng(), rng());
                stim::simd_word rx2(rng(), rng());
                stim::simd_word rz2(rng(), rng());

                // Only update phase if the qubits are not leaked.
                r ^= (x1 & z2) & ~(x2 ^ z1) & ~(l1 | l2);
                x2 ^= ~(l1 | l2) & x1;
                z1 ^= ~(l1 | l2) & z2;
                // Apply depolarizing error if they are leaked.
                r ^= (z1 & rx1 & l2)
                        & (x1 & rz1 & l2)
                        & (z2 & rx2 & l1)
                        & (x2 & rz2 & l2);
            });
        }
    }
}

void
CliffordSimulator::M(std::vector<uint> operands, bool record) {
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
            // Predicate the rowsum where xij = 1
            uint k = n_qubits*i + j;
            browsum(2*n_qubits, i+n_qubits, true, x_table[k]);
        }
        // Measurement outcome is r_table[2*n_qubits]
        // Condition results based on x_table[x_width-1] = 0
        // Also, if any of these qubits were leaked, set the measurement
        // output to a random value.
        r_table[2*n_qubits].for_each_word(
                x_table[x_width-1], leak_table[j],
        [&] (auto& r, auto& x, auto& l)
        {
            stim::simd_word rand(rng(), rng());
            r = (r & ~x & ~l) | (rand & l);
        });
        // Now, perform any indeterminate measurements.
        // This takes more time as we must handle each case separately.
        for (uint64_t t = 0; t < shots; t++) {
            if (!x_table[x_width-1][t]) continue;
            // Find first i where xia = 1 in i = n+1 to 2n. Call this ii.
            uint ii;
            for (uint i = n_qubits; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (x_table[k][t]) {
                    ii = i;
                    break;
                }
            }
            // First perform rowsums.
            for (uint i = 0; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (i != ii && x_table[k][t]) {
                    rowsum(i, ii, t);
                }
            }
            // Second, swap (ii-n)-th row with ii-th row and clear
            // ii-th row. Set rii to random value.
            for (uint jj = 0; jj < n_qubits; jj++) {
                x_table[(ii-n_qubits)*n_qubits+jj][t] = x_table[ii*n_qubits+jj][t];
                z_table[(ii-n_qubits)*n_qubits+jj][t] = z_table[ii*n_qubits+jj][t];
                x_table[ii*n_qubits+jj][t] = 0;
                z_table[ii*n_qubits+jj][t] = 0;
            }
            r_table[ii-n_qubits][t] = r_table[ii][t];
            r_table[ii][t] = rng() & 1;
            z_table[ii*n_qubits+j][t] = 1;
            // For simplicity, copy the rii value (the measurement outcome)
            // to the scratch space (row 2n+1) in r.
            r_table[2*n_qubits][t] = r_table[ii][t];
        }

        if (record) {
            record_table[record_offset].clear();
            record_table[record_offset++].swap_with(r_table[2*n_qubits]);
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
    uint64_t r = record_offset;
    M(operands, true);
    while (record_offset > r) {
        uint j = operands[record_offset-r-1];
        for (uint i = 0; i < 2*n_qubits; i++) {
            uint k = i*n_qubits + j;
            r_table[i].for_each_word(
                z_table[k], record_table[record_offset-1],
            [&](auto& r, auto& z, auto& rec) {
                r ^= z & rec;
            });
        }
        record_table[--record_offset].clear();
    }
}

void
CliffordSimulator::eDPO(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            for (uint ii = 0; ii < 2*n_qubits; ii++) {
                uint k = n_qubits*ii + j;
                // Implement this as a reset.
                x_table[k][t] &= z_table[k][t];
                r_table[ii][t] ^= x_table[k][t];
                x_table[k][t] = 0;
                z_table[k][t];
                // Not sure if |2> also relaxes to |0> at the same rate. 
                // We'll make a conservative assumption that if the qubit
                // is leaked, it is now in |1>
                r_table[ii][t] ^= leak_table[j][t];
                leak_table[j][t] = 0;
            }
        });
    }
}

void
CliffordSimulator::eDPH(std::vector<uint> operands, std::vector<fp_t> rates) {
    // I believe that properly simulating dephasing errors requires a state
    // vector simulator, so to a first order approximation, we will just
    // perform a "Hadamard error" so we rotate the qubit by 90 deg. This
    // results in a uniform superposition in the current basis.
    //  i.e. |+> --> |0> = |+> + |->    (ignoring the constant for readability)
    //       |0> --> |+> = |0> + |1>
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            for (uint ii = 0; ii < 2*n_qubits; ii++) {
                uint k = n_qubits*ii + j;
                r_table[ii][t] ^= x_table[k][t] & z_table[k][t] & ~leak_table[j][t];
                auto tmp = z_table[k][t];
                z_table[k][t] = (z_table[k][t] & leak_table[j][t])
                                | (x_table[k][t] & ~leak_table[j][t]);
                x_table[k][t] = (x_table[k][t] & leak_table[j][t])
                                | (tmp & ~leak_table[j][t]);
            }
        });
    }
}

void
CliffordSimulator::eDP1(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            auto p = rng() & 3;
            for (uint ii = 0; ii < 2*n_qubits; ii++) {
                uint k = n_qubits*ii + j;
                r_table[ii][t] ^= ((z_table[k][t]) & (p & 1))
                                    & ((x_table[k][t]) & (p & 2))
                                    & ~leak_table[j][t];
            }
        });
    }
}

void
CliffordSimulator::eDP2(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i];
        uint j2 = operands[i+1];

        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) {
            auto p = rng() & 15;
            for (uint ii = 0; ii < 2*n_qubits; ii++) {
                uint k1 = n_qubits*ii + j1;
                uint k2 = n_qubits*ii + j2;
                r_table[ii][t] ^= ((z_table[k1][t]) & (p & 1))
                                    & ((x_table[k1][t]) & (p & 2))
                                    & ~leak_table[j1][t];
                r_table[ii][t] ^= ((z_table[k2][t]) & (p & 4))
                                    & ((x_table[k2][t]) & (p & 8))
                                    & ~leak_table[j2][t];
                // Todo: flag correlated errors
            }
        });
    }
}

void
CliffordSimulator::eX(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            for (uint ii = 0; ii < 2*n_qubits; ii++) {
                uint k = n_qubits*ii + j;
                r_table[ii][t] ^= z_table[k][t] & ~leak_table[j][t];
            }
        });
    }
}

void
CliffordSimulator::eLI(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i];
        uint j2 = operands[i+1];

        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) {
            auto p = rng() % 3;
            // Either leak (or unleak) one or both qubits.
            leak_table[j1][t] ^= (p == 0) || (p == 2);
            leak_table[j2][t] ^= (p == 1) || (p == 2);
        });
    }
}

void
CliffordSimulator::eLT(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i];
        uint j2 = operands[i+1];

        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) {
            leak_table[j1][t] |= leak_table[j2][t];
            leak_table[j2][t] |= leak_table[j1][t];
        });
    }
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
    x_table_cpy = stim::simd_bit_table(x_table);
    z_table_cpy = stim::simd_bit_table(z_table);
    r_table_cpy = stim::simd_bit_table(r_table);
    leak_table_cpy = stim::simd_bit_table(leak_table);
}

void
CliffordSimulator::rollback_at_trial(uint64_t t) {
    for (uint i = 0; i < x_width; i++)  x_table[i][t] = x_table_cpy[i][t];
    for (uint i = 0; i < z_width; i++)  z_table[i][t] = z_table_cpy[i][t];
    for (uint i = 0; i < r_width; i++)  r_table[i][t] = r_table_cpy[i][t];
    for (uint i = 0; i < leak_width; i++)  leak_table[i][t] = leak_table_cpy[i][t];
}

}   // qontra
