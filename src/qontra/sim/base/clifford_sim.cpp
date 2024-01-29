/*
 *  author: Suhas Vittal
 *  date:   18 June 2023
 * */

#include "qontra/sim/base/clifford_sim.h"
#include "qontra/sim/base/gate_impl/clifford_sim.h"

namespace qontra {

inline uint64_t csim_x_width(uint64_t n) { return 2*n*(n+1) + 1; }
inline uint64_t csim_z_width(uint64_t n) { return 2*n*(n+1); }
inline uint64_t csim_r_width(uint64_t n) { return 2*n + 4; }

CliffordSimulator::CliffordSimulator(uint64_t n, uint64_t max_shots)
    :StateSimulator(n, max_shots),
    x_table(csim_x_width(n), max_shots),
    z_table(csim_z_width(n), max_shots),
    r_table(csim_r_width(n), max_shots),
    x_table_cpy(csim_x_width(n), max_shots),
    z_table_cpy(csim_z_width(n), max_shots),
    r_table_cpy(csim_r_width(n), max_shots),
    x_width(csim_x_width(n)),
    z_width(csim_z_width(n)),
    r_width(csim_r_width(n))
{
    init_tables();
}

CliffordSimulator::CliffordSimulator(const CliffordSimulator& other)
    :StateSimulator(other),
    x_table(other.x_table),
    z_table(other.z_table),
    r_table(other.r_table),
    x_table_cpy(other.x_table_cpy),
    z_table_cpy(other.z_table_cpy),
    r_table_cpy(other.r_table_cpy),
    x_width(other.x_width),
    z_width(other.z_width),
    r_width(other.r_width)
{}

void
CliffordSimulator::H(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (uint64_t q : operands) {
            uint64_t k = get_index(i, q);
            if (tr >= 0) {
                stim::bit_ref r = r_table[i][tr],
                                x = x_table[k][tr],
                                z = z_table[k][tr],
                                lock = lock_table[q][tr];
                __h_gate(r, x, z, lock);
            } else {
                r_table[i].for_each_word(
                    x_table[k], z_table[k], lock_table[q],
                [](auto& r, auto& x, auto& z, auto& lock)
                {
                    __h_gate(r, x, z, lock);
                });
            }
        }
    }
}

void
CliffordSimulator::X(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (uint64_t q : operands) {
            uint64_t k = get_index(i, q);
            if (tr >= 0) {
                stim::bit_ref r = r_table[i][tr],
                                z = z_table[k][tr],
                                lock = lock_table[q][tr];
                __x_gate(r, z, lock);
            } else {
                r_table[i].for_each_word(
                    z_table[k], lock_table[q],
                [&](auto& r, auto& z, auto& lock)
                {
                    __x_gate(r, z, lock);
                });
            }
        }
    }
}

void
CliffordSimulator::Z(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (uint64_t q : operands) {
            uint64_t k = get_index(i, q);
            if (tr >= 0) {
                stim::bit_ref r = r_table[i][tr],
                                x = x_table[k][tr],
                                lock = lock_table[q][tr];
                __z_gate(r, x, lock);
            } else {
                r_table[i].for_each_word(
                    x_table[k], lock_table[q],
                [&](auto& r, auto& x, auto& lock)
                {
                    __z_gate(r, x, lock);
                });
            }
        }
    }
}

void
CliffordSimulator::S(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (uint64_t q : operands) {
            uint64_t k = get_index(i, q);
            if (tr >= 0) {
                stim::bit_ref r = r_table[i][tr],
                                x = x_table[k][tr],
                                z = z_table[k][tr],
                                lock = lock_table[q][tr];
                __s_gate(r, x, z, lock);
            } else {
                r_table[i].for_each_word(
                        x_table[k], z_table[k], lock_table[q],
                [&](auto& r, auto& x, auto& z, auto& lock)
                {
                    __s_gate(r, x, z, lock);
                });
            }
        }
    }
}

void
CliffordSimulator::CX(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (size_t j = 0; j < operands.size(); j += 2) {
            uint64_t q1 = operands[j],
                     q2 = operands[j+1];
            uint64_t k1 = get_index(i, q1),
                     k2 = get_index(i, q2);

            stim::simd_bits lock_both(lock_table[q1]);
            lock_both |= lock_table[q2];

            if (tr >= 0) {
                stim::bit_ref r = r_table[i][tr],
                                x1 = x_table[k1][tr],
                                z1 = z_table[k1][tr],
                                x2 = x_table[k2][tr],
                                z2 = z_table[k2][tr],
                                lock = lock_both[tr];
                __cx_gate(r, x1, z1, x2, z2, lock);
            } else {
                std::array<stim::simd_bits_range_ref<SIMD_WIDTH>, 6> args{{
                    r_table[i],
                    x_table[k1],
                    z_table[k1],
                    x_table[k2],
                    z_table[k2],
                    lock_both
                }};
                for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>*, 6>& word_array)
                {
                    auto& r = *word_array[0],
                        & x1 = *word_array[1],
                        & z1 = *word_array[2],
                        & x2 = *word_array[3],
                        & z2 = *word_array[4],
                        & lock = *word_array[5];
                    __cx_gate(r, x1, z1, x2, z2, lock);
                });
            }
        }
    }
}

void
CliffordSimulator::M(
        std::vector<uint64_t> operands,
        std::vector<fp_t> m1w0,
        std::vector<fp_t> m0w1,
        int record,
        int64_t tr)
{
    size_t k = 0;
    for (uint64_t q : operands) {
        clear_scratch_space();
        if (tr >= 0) {
            measure_qubit_in_trial(q, tr);
            if (record >= 0 && !lock_table[q][tr]) {
                record_table[record][tr] = r_table[2*n_qubits][tr];
                fp_t e = get_probability_sample_from_rng();
                if (record_table[record][tr] == 0) {
                    if (e < m1w0[k]) record_table[record][tr] = 1;
                } else {
                    if (e < m0w1[k]) record_table[record][tr] = 1;
                }
                record++;
                k++;
            }
        } else {
            batch_measure_qubit(q);
            if (record >= 0) {
                record_table[record].for_each_word(r_table[2*n_qubits], lock_table[q],
                        [&] (auto& rec, auto& r, auto& lock)
                        {
                            rec |= andnot(lock, r);
                        });
                r_table[2*n_qubits].clear();
                error_channel_m(record, m1w0[k], m0w1[k], lock_table[q]);
                record++;
                k++;
            }
        }
    }
}

void
CliffordSimulator::R(std::vector<uint64_t> operands, int64_t tr) {
    // Implement as a measure + X gate.
    stim::simd_bit_table<SIMD_WIDTH> record_cpy = record_table;
    std::vector<fp_t> all_zeros(operands.size(), 0.0);

    M(operands, all_zeros, all_zeros, 0, tr);
    for (uint64_t i = 0; i < 2*n_qubits; i++) {
        for (size_t j = 0; j < operands.size(); j++) {
            uint64_t q = operands[j];
            uint64_t k = get_index(i, q);

            if (tr >= 0) {
                r_table[i][tr] &= z_table[k][tr] & record_table[j][tr] & lock_table[q][tr];
            } else {
                r_table[i].for_each_word(z_table[k], record_table[j], lock_table[q],
                        [] (auto& r, auto& z, auto& rec, auto& lock)
                        {
                            r ^= z & rec & lock;
                        });
            }
        }
    }
    // Recover the record table entries that we had used.
    for (size_t j = 0; j < operands.size(); j++) {
        if (tr >= 0) {
            record_table[j][tr].swap_with(record_cpy[j][tr]);
        } else {
            record_table[j].swap_with(record_cpy[j]);
        }
    }
}

void
CliffordSimulator::rowsum(uint64_t h, uint64_t i, uint64_t t) {
    // Clear scratch storage:
    r_table[2*n_qubits+1][t] = 0;
    r_table[2*n_qubits+2][t] = 0;

    r_table[2*n_qubits+2][t] = r_table[h][t] ^ r_table[i][t];
    for (size_t j = 0; j < n_qubits; j++) {
        uint64_t kh = get_index(h, j),
                 ki = get_index(i, j);

        stim::bit_ref x1 = x_table[ki][t],
                        x2 = x_table[kh][t],
                        z1 = z_table[ki][t],
                        z2 = z_table[kh][t],
                        s1 = r_table[2*n_qubits+1][t],
                        s2 = r_table[2*n_qubits+2][t];
        __rowsum_arith(x1, z1, x2, z2, s1, s2, true);
    }
    r_table[h][t] = r_table[2*n_qubits+2][t];
}

void
CliffordSimulator::browsum(uint64_t h, uint64_t i, stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    // Clear the scratch storage used.
    r_table[2*n_qubits+1].clear();
    r_table[2*n_qubits+2].clear();
    // 2n + 1 and 2n + 2 correspond to the 1s and 2s place in a 
    // 2-bit number.
    //
    // To perform the 2rh + 2ri computation, we can just
    // XOR the bits going to the 2s place (as this is mod 4).
    r_table[2*n_qubits+2] ^= r_table[h];
    r_table[2*n_qubits+2] ^= r_table[i];
    for (uint64_t j = 0; j < n_qubits; j++) {
        uint64_t kh = get_index(h, j),
                 ki = get_index(i, j);
        std::array<stim::simd_bits_range_ref<SIMD_WIDTH>, 7> args{{
            x_table[ki],
            z_table[ki],
            x_table[kh],
            z_table[kh],
            r_table[2*n_qubits+1],
            r_table[2*n_qubits+2],
            pred
        }};
        for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>*, 7>& word_array)
        {
            auto& x1 = *word_array[0],
                & z1 = *word_array[1],
                & x2 = *word_array[2],
                & z2 = *word_array[3],
                & s1 = *word_array[4],
                & s2 = *word_array[5],
                & p = *word_array[6];
            __rowsum_arith(x1, z1, x2, z2, s1, s2, p);
        });
    }   // After all the iterations, row 2n+1 and 2n+2 should have 
        // the value of the rowsum.

    // Now, in locations where row 2n+2 == 0, we need to set rh = 0. Otherwise,
    // we set it to 1. This is as simple as just moving row 2n+2 to rh.
    r_table[h].for_each_word(r_table[2*n_qubits+2], pred,
    [&] (auto& rh, auto& s, auto& p) {
        rh = (p & s) | andnot(p, rh);
    });
}

void
CliffordSimulator::batch_measure_qubit(uint64_t q) {
    // We need to check amongst the stabilizer rows (n+1 to 2n) whether there is a p such that
    // x_pq = 1.
    //
    // This is easily done by simply ORing all rows in the range.
    for (uint64_t i = n_qubits; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        x_table[x_width-1] |= x_table[k];
    }
    // If the qubit is locked, treat the measurement as deterministic.
    lock_table[q].invert_bits();
    x_table[x_width-1] &= lock_table[q];
    // As the deterministic case is easy to handle, we will do this first (where x_pq = 0).
    // 
    // First, clear row 2*n+1.
    clear_row(2*n_qubits);
    // Perform rowsum between row 2*n+1 and all rows i+n for 1 <= i <= n where x_iq = 1.
    for (uint64_t i = 0; i < n_qubits; i++) {
        uint64_t k = get_index(i, q);
        // Predicate this rowsum where x_ij = 1 and j is unlocked.
        stim::simd_bits<SIMD_WIDTH> pred(x_table[k]);
        pred &= lock_table[q];
        browsum(2*n_qubits, i+n_qubits, pred);
    }
    lock_table[q].invert_bits();
    // The measurement results are in r_table[2*n_qubits].
    //
    // Filter the measurement results to meet the x_pj = 0 condition.
    r_table[2*n_qubits].for_each_word(x_table[x_width-1], lock_table[q],
        [] (auto& r, auto& x, auto& lock) {
            r = andnot(x | lock, r);
        });
    // Now, handle the indeterminate case.
    //
    // We will try to batch this as much as possible.
    stim::simd_bits completed_trials(x_table[x_width-1]);
    completed_trials.invert_bits(); // Now, = 1 wherever a measurement is already done.
    for (uint64_t t = 0; t < shots; t++) {
        if (completed_trials[t]) continue;
        // Find first p where x_pq = 1 in p = n+1 to 2n.
        uint64_t p;
        for (uint64_t i = n_qubits; i < 2*n_qubits; i++) {
            uint64_t k = get_index(i, q);
            if (x_table[k][t]) { p = i; break; }
        }
        // Collect all trials where x_pq = 1.
        stim::simd_bits<SIMD_WIDTH> trials_with_good_p(x_table[get_index(p, q)]);
        // Set all completed trials in trials_with_good_p to 0.
        completed_trials.invert_bits();
        trials_with_good_p &= completed_trials;
        completed_trials.invert_bits();
        // Now, perform the indeterminate measurement.
        //
        // First, do a rowsum on all rows i where x_iq = 1 and i != p.
        for (uint64_t i = 0; i < 2*n_qubits; i++) {
            uint64_t k = get_index(i, q);
            if (i == p) continue;
            // Create browsum predicate.
            stim::simd_bits<SIMD_WIDTH> pred(trials_with_good_p);
            pred &= x_table[k];
            browsum(i, p, pred);
        }
        // Now, swap row p-n with row p and clear the p-th row. Set r_p to a random value.
        //
        // Also, for simplicity, copy r_p to the 2n+1 row (2n+1 row stores measurement outcomes for
        // the qubit).
        swap_rows_where(p-n_qubits, p, trials_with_good_p);
        clear_row_where(p, trials_with_good_p);
        stim::simd_bits<SIMD_WIDTH> rand_meas =
            stim::simd_bits<SIMD_WIDTH>::random(r_table.num_minor_bits_padded(), rng);
        r_table[p].for_each_word(r_table[2*n_qubits], rand_meas, trials_with_good_p,
                [] (auto& r1, auto& r2, auto& rand, auto& tp)
                {
                    r1 = (tp & rand) | andnot(tp, r1);
                    r2 = (tp & r1) | andnot(tp, r2);
                });
        z_table[get_index(p, q)] |= trials_with_good_p;
        // Update completed_trials.
        completed_trials |= trials_with_good_p;
        if (completed_trials.popcnt() == shots) break;
    }
}

void
CliffordSimulator::measure_qubit_in_trial(uint64_t q, int64_t tr) {
    if (lock_table[q][tr]) return;  // Cannot measure the qubit.
    // Simply perform the algorithm in batch_measure_qubit, but optimize it for a single trial.
    for (uint64_t i = n_qubits; i < 2*n_qubits; i++) {
        uint64_t k = get_index(i, q);
        x_table[x_width-1][tr] |= x_table[k][tr];
    }
    if (x_table[x_width-1][tr]) {
        uint64_t p;
        for (uint64_t i = n_qubits; i < 2*n_qubits; i++) {
            uint64_t k = get_index(i, q);
            if (x_table[k][tr]) { p = i; break; }
        }
        for (uint64_t i = 0; i < 2*n_qubits; i++) {
            uint64_t k = get_index(i, q);
            if (i != p && x_table[k][tr]) rowsum(i, p, tr);
        }
        swap_rows_for_trial(p-n_qubits, p, tr);
        clear_row_for_trial(p, tr);
        r_table[p][tr] = rng() & 1;
        z_table[get_index(p, q)][tr] = 1;
        r_table[2*n_qubits][tr] = r_table[p][tr];
    } else {
        // Complete deterministic measurement.
        clear_row_for_trial(2*n_qubits, tr);
        for (uint64_t i = 0; i < n_qubits; i++) {
            uint64_t k = get_index(i, q);
            if (x_table[k][tr]) rowsum(2*n_qubits, i+n_qubits, tr);
        }
    }
}

}   // qontra
