/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#include "sim.h"

namespace contra {

// Default values
#define SIM_MAX_QUBITS  1024
#define SIM_T(p)       (-1000.0/log(1-(p)) * 1e3)

ContraSim::ContraSim(uint distance, fp_t p, Code qec_code)
    :error_params(SIM_MAX_QUBITS, p),
    time_params(SIM_MAX_QUBITS, 30, 40, SIM_T(p), SIM_T(p)),
    distance(distance),
    qec_code(qec_code),
    qc_x_table(1, 1),
    qc_z_table(1, 1),
    qc_syndrome_buf(1, 1),
    cp_syndrome_history(1, 1),
    cp_pauli_frames(1, 1),
    cp_pauli_frames_ideal(1, 1),
    pc(0),
    register_file(),
    imem(),
    dmem()
{
    register_file.fill((Register) {0, false});
    dmem.fill(0);
}

void
ContraSim::qc_blk_recv_cmd() {
    // Edit this region if you need to modify syndrome extraction.
    for (auto& inst : instruction_buf) {
        stim::simd_bit_table x_cpy(qc_x_table_rwidth, 
                inst.exclude_trials.size());
        stim::simd_bit_table z_cpy(qc_z_table_rwidth, 
                inst.exclude_trials.size());
        stim::simd_bit_table r_cpy(qc_r_table_rwidth, 
                inst.exclude_trials.size());
        stim::simd_bit_table leak_cpy(qc_leak_table_rwidth, 
                inst.exclude_trials.size());

        // Copy X and Z tables
        for (uint j = 0; j < xz_size; j++) {
            for (uint i = 0; i < inst.exclude_trials.size(); i++) {
                uint64_t t = inst.exclude_trials[i];
                x_cpy[j][t] ^= qc_x_table[j][t];
                z_cpy[j][t] ^= qc_z_table[j][t];
            }
        }
        // Copy R table
        for (uint j = 0; j < r_size; j++) {
            for (uint i = 0; i < inst.exclude_trials.size(); i++) {
                uint64_t t = inst.exclude_trials[i];
                r_cpy[j][t] ^= qc_r_table[j][t];
            }
        }
        // Copy leakage table
        for (uint j = 0; j < n_qubits; j++) {
            for (uint i = 0; i < inst.exclude_trials.size(); i++) {
                uint64_t t = inst.exclude_trials[i];
                leak_cpy[j][t] ^= qc_leak_table[j][t];
            }
        }

        switch (inst.name) {
        case "H":
            qc_H(inst.operands);
            break;
        case "X":
            qc_X(inst.operands);
            break;
        case "Z":
            qc_Z(inst.operands);
            break;
        case "S":
            qc_S(inst.operands);
            break;
        case "CX":
            qc_CX(inst.operands);
            break;
        case "Mnrc":
            qc_M(inst.operands, false);
            break;
        case "Mrc":
            qc_M(inst.operands, true);
            break;
        case "R":
            qc_R(inst.operands);
            break;
        }
    }
}

void
ContraSim::qc_H(const std::vector<uint>& operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            qc_r_table[i].for_each_word(
                qc_x_table[k], qc_z_table[k], qc_leak_table[j],
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
ContraSim::qc_X(const std::vector<uint>& operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            // Only flip where not leaked.
            qc_leak_table[j].invert_bits();  // now non-leaked locations are 1.
            qc_x_table[k] ^= qc_leak_table[j];
            qc_leak_table[j].invert_bits();
        }
    }
}

void
ContraSim::qc_Z(const std::vector<uint>& operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            // Only flip where not leaked.
            qc_leak_table[j].invert_bits();  // now non-leaked locations are 1.
            qc_z_table[k] ^= qc_leak_table[j];
            qc_leak_table[j].invert_bits();
        }
    }
}

void
ContraSim::qc_S(const std::vector<uint>& operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = i*n_qubits + j;
            qc_r_table[i].for_each_word(
                    qc_x_table[k], qc_z_table[k], qc_leak_table[j],
            [&](auto& r, auto& x, auto& z, auto& l)
            {
                r ^= x & z & ~l;
                z ^= x & ~l;
            });
        }
    }
}

void
ContraSim::qc_CX(const std::vector<uint>& operands) {
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint ii = 0; ii < operands.size(); ii += 2) {
            uint j1 = inst.operands[ii];
            uint j2 = inst.operands[ii+1];
            uint k1 = i*n_qubits + j1;
            uint k2 = i*n_qubits + j2;

            qc_r_table[i].for_each_word(
                qc_x_table[k1], qc_z_table[k1],
                qc_x_table[k2], qc_z_table[k2],
                qc_leak_table[j1], qc_leak_table[j2],
            [&](auto& r, auto& x1, auto& z1, auto& x2, auto& z2, auto& l1, auto& l2)
            {
                stim::simd_word rx1(rng(), rng());
                stim::simd_word rz1(rng(), rng());
                stim::simd_word rx2(rng(), rng());
                stim::simd_word rz2(rng(), rng());

                // Only update phase if the qubits are not leaked.
                // Otherwise, apply a random depolarizing error.
                r ^= (x1 & z2) & ~(x2 ^ z1) & ~(l1 | l2);

                x2 ^= (l1 & rx2) | (~l1 & x1);
                z1 ^= (l2 & rz1) | (~l2 & z1);
                x1 ^= l2 & rx1;
                z2 ^= l1 & rz2;
            });
        }
    }
}

void
ContraSim::qc_M(const std::vector<uint>& operands, bool record_in_syndrome_buf) {
    for (uint j : operands) {
        // Clear scratch space
        qc_x_table[2*n_qubits].clear();
        qc_r_table[2*n_qubits].clear();
        qc_r_table[2*n_qubits+1].clear();
        qc_r_table[2*n_qubits+2].clear();
        qc_r_table[2*n_qubits+3].clear();

        // Figure out if xij is nonzero anywhere (just need to OR all rows).
        for (uint i = n_qubits+1; i < 2*n_qubits; i++) {
            uint k = n_qubits*i + j;
            qc_x_table[2*n_qubits] |= qc_x_table[k];
        }
        // Also exclude rows where the qubit is leaked.
        qc_leak_table[j].invert_bits();
        qc_x_table[2*n_qubits] &= qc_leak_table[j];
        qc_leak_table[j].invert_bits();
        
        // Perform determinate measurement. We want to condition this
        // to happen wherever qc_x_table[2*n_qubits] = 1
        for (uint i = 0; i < n_qubits; i++) {
            // Predicate the rowsum where xij = 1
            uint k = n_qubits*i + j;
            rowsum(2*n_qubits, i+n_qubits, true, qc_x_table[k]);
        }
        // Measurement outcome is qc_r_table[2*n_qubits]
        // Condition results based on qc_x_table[2*n_qubits] = 0
        // Also, if any of these qubits were leaked, set the measurement
        // output to a random value.
        qc_r_table[2*n_qubits].for_each_word(
                qc_x_table[2*n_qubits], qc_leak_table[j],
        [&] (auto& r, auto& x, auto& l)
        {
            stim::simd_word rand(rng(), rng());
            r = (r & ~x & ~l) | (rand & l);
        });
        // Now, perform any indeterminate measurements.
        // This takes more time as we must handle each case separately.
        for (uint64_t t = 0; t < n_trials; t++) {
            if (!qc_x_table[2*n_qubits][t]) continue;
            // Find first i where xia = 1 in i = n+1 to 2n. Call this ii.
            uint ii;
            for (uint i = n_qubits+1; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (qc_x_table[k][t]) {
                    ii = i;
                    break;
                }
            }
            // First perform rowsums.
            for (uint i = 0; i < 2*n_qubits; i++) {
                uint k = n_qubits*i + j;
                if (i != ii && qc_x_table[k][t]) {
                    rowsum(i, ii, false, qc_x_table[0]);
                }
            }
            // Second, swap (ii-n)-th row with ii-th row and clear
            // ii-th row. Set rii to random value.
            qc_r_table[ii-n_qubits][t] = qc_r_table[ii][t];
            qc_r_table[ii][t] = rng() & 1;
            for (uint jj = 0; jj < n_qubits; jj++) {
                uint k1 = n_qubits*(ii-n_qubits) + jj;
                uint k2 = n_qubits*ii + jj;
                qc_x_table[k1][t] = qc_x_table[k2][t];
                qc_x_table[k2][t] = 0;
            }
            // For simplicity, copy the rii value (the measurement outcome)
            // to the scratch space (row 2n+1) in r.
            qc_r_table[2*n_qubits][t] = qc_r_table[ii][t];
        }
        // Drop the measurement in the syndrome buffer
        if (record_in_syndrome_buf) {
            qc_syndrome_buf[syndrome_buf_offset] |= qc_r_table[2*n_qubits];
    `       if (syndrome_buf_lookback > 0) {
                qc_syndrome_buf[syndrome_buf_offset] ^=
                    qc_syndrome_buf[syndrome_buf_offset - syndrome_buf_lookback];
            }
            syndrome_buf_offset++;
            next_buf_lookback++;
        }
    }
}

void
ContraSim::qc_R(const std::vector<uint>& operands) {
    // We can just implement this as a "CNOT" between a qubit and itself.
    for (uint i = 0; i < 2*n_qubits; i++) {
        for (uint j : operands) {
            uint k = n_qubits*i + j;
            qc_x_table[k] &= qc_z_table[k];
            qc_r_table[i] ^= qc_x_table[k];
            qc_x_table[k].clear();
            qc_z_table[k].clear();
        }
    }
}

void
ContraSim::rowsum(uint h, uint i, bool use_pred, stim::simd_bit_range_ref& pred) {
    // 2n + 1 and 2n + 2 correspond to the 1s and 2s place in a 
    // 2-bit number.
    //
    // To perform the 2rh + 2ri computation, we can just
    // XOR the bits going to the 2s place (as this is mod 4).
    const stim::simd_word en_p(0, 0);
    if (use_pred)   en_p = ~en_p;

    qc_r_table[2*n_qubits+2] ^= qc_r_table[h] ^ qc_r_table[i];
    for (uint j = 0; j < n_qubits; j++) {
        uint kh = n_qubits * h + j;
        uint ki = n_qubits * i + j;
        qc_x_table[ki].for_each_word(
            qc_z_table[ki], qc_x_table[kh], qc_z_table[kh],
            qc_r_table[2*n_qubits+1], qc_r_table[2*n_qubits+2],
            pred
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
            stim::simd_word gsum_mag(0, 0);
            stim::simd_word gsum_sgn(0, 0);jkkkkj
            
            gsum_mag ^= (~x1 & z1 & x2)
                        | (x1 & ~z1 & z2) | (x1 & z1 & (x2 ^ z2));
            gsum_sgn ^= (~x1 & z1 & x2 & x2) 
                        | (x1 & ~z1 & x2 & z2)
                        | (x1 & z1 & x2 & ~z2);
            // Perform the modular addition and save the results in
            // the scratch storage
            scr2 ^= (en_p & p) & ((scr1 & gsum_mag) ^ gsum_sgn);
            scr1 ^= (en_p & p) & gsum_mag;
            // Update x2 and z2
            x2 ^= (en_p & p) & x1;
            z2 ^= (en_p & p) & z1;
        }
    }
    // Now, in locations where row 2n+2 == 0, we need to set rh = 0. Otherwise,
    // we set it to 1. This is as simple as just moving row 2n+2 to rh.
    qc_r_table[h].for_each_word(qc_r_table[2*n_qubits+2], pred
    [&] (auto& rh, auto& scr, auto& p) {
        rh = ((en_p & p) & scr) | (~(en_p & p) & rh);
    });
    // Clear the scratch storage used.
    qc_r_table[2*n_qubits+1].clear();
    qc_r_table[2*n_qubits+2].clear();
}

}   // contra

