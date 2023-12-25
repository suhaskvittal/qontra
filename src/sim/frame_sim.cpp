/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "sim/frame_sim.h"
#include "sim/frame_sim_gate_impl.h"

namespace qontra {

#define L_I     StateSimulator::label_1q_t::I
#define L_X     StateSimulator::label_1q_t::X
#define L_Y     StateSimulator::label_1q_t::Y
#define L_Z     StateSimulator::label_1q_t::Z
#define L_L     StateSimulator::label_1q_t::L

static const StateSimulator::label_1q_t IXYZ[4] = { L_I, L_X, L_Y, L_Z };
static const StateSimulator::label_1q_t IL[2] = { L_I, L_L };

void
FrameSimulator::H(std::vector<uint> operands, int64_t tr) {
    for (uint i : operands) {
        if (tr >= 0) {
            stim::bit_ref x = x_table[i][tr],
                            z = z_table[i][tr],
                            lock = lock_table[i][tr];
            __h_gate(x, z, lock);
        } else {
            x_table[i].for_each_word(
                    z_table[i],
                    lock_table[i],
                    [] (auto& x, auto& z, auto& lock)
                    {
                        __h_gate(x, z, lock);
                    });
        }
    }
}

void
FrameSimulator::X(std::vector<uint> operands, int64_t tr) {
    for (uint i : operands) {
        // Flip X wherever lock = 0
        if (tr >= 0) {
            stim::bit_ref x = x_table[i][tr],
                            lock = lock_table[i][tr];
            __x_gate(x, lock);
        } else {
            lock_table[i].invert_bits();
            x_table[i] ^= lock_table[i];
            lock_table[i].invert_bits();
        }
    }
}

void
FrameSimulator::Z(std::vector<uint> operands, int64_t tr) {
    for (uint i : operands) {
        // Flip Z wherever lock = 0
        if (tr >= 0) {
            stim::bit_ref z = z_table[i][tr],
                            lock = lock_table[i][tr];
            __z_gate(z, lock);
        } else {
            lock_table[i].invert_bits();
            z_table[i] ^= lock_table[i];
            lock_table[i].invert_bits();
        }
    }
}

void
FrameSimulator::S(std::vector<uint> operands, int64_t tr) {} // Do not implement.

void
FrameSimulator::CX(std::vector<uint> operands, int64_t tr) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
    
        stim::simd_bits<SIMD_WIDTH> lock_both(lock_table[j1]);
        lock_both |= lock_table[j2];

        if (tr >= 0) {
            uint64_t rand_word = rng();
            stim::bit_ref x1 = x_table[j1][tr],
                            z1 = z_table[j1][tr],
                            l1 = leak_table[j1][tr],
                            x2 = x_table[j2][tr],
                            z2 = z_table[j2][tr],
                            l2 = leak_table[j2][tr],
                            lock = lock_both[tr];
            bool rx1 = rand_word & 0x1,
                 rz1 = rand_word & 0x2,
                 rx2 = rand_word & 0x4,
                 rz2 = rand_word & 0x8;
            __cx_gate(x1, z1, l1, x2, z2, l2, lock, rx1, rz1, rx2, rz2);
        } else {
            const size_t width = x_table[j1].num_bits_padded();
            stim::simd_bits<SIMD_WIDTH> rand_x_1 = stim::simd_bits<SIMD_WIDTH>::random(width, rng),
                                        rand_z_1 = stim::simd_bits<SIMD_WIDTH>::random(width, rng),
                                        rand_x_2 = stim::simd_bits<SIMD_WIDTH>::random(width, rng),
                                        rand_z_2 = stim::simd_bits<SIMD_WIDTH>::random(width, rng);
            std::array<stim::simd_bits_range_ref<SIMD_WIDTH>, 11> args{{
                x_table[j1],
                z_table[j1],
                leak_table[j1],
                x_table[j2],
                z_table[j2],
                leak_table[j2],
                lock_both,
                rand_x_1,
                rand_z_1,
                rand_x_2,
                rand_z_2
            }};
            for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>, 11> word_array) {
                auto& x1 = word_array[0],
                    & z1 = word_array[1],
                    & l1 = word_array[2],
                    & x2 = word_array[3],
                    & z2 = word_array[4],
                    & l2 = word_array[5],
                    & lock = word_array[6],
                    & rx1 = word_array[7],
                    & rz1 = word_array[8],
                    & rx2 = word_array[9],
                    & rz2 = word_array[10];
                __cx_gate(x1, z1, l1, x2, z2, l2, lock, rx1, rz1, rx2, rz2);
            });
        }
    }
}

void
FrameSimulator::LEAKAGE_ISWAP(std::vector<uint> operands, int64_t tr) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
    
        stim::simd_bits<SIMD_WIDTH> lock_both(lock_table[j1]);
        lock_both |= lock_table[j2];

        if (tr >= 0) {
            stim::bit_ref x1 = x_table[j1][tr],
                            x2 = x_table[j2][tr],
                            l2 = leak_table[j2][tr],
                            lock = lock_both[tr];
            __liswap_gate(x1, x2, l2, lock);
        } else {
            x_table[j1].for_each_word(
                x_table[j2],
                leak_table[j2],
                lock_both,
            [] (auto& x1, auto& x2, auto& l2, auto& lock)
            {
                __liswap_gate(x1, x2, l2, lock);
            });
        }
    }
}

void
FrameSimulator::M(
        std::vector<uint> operands, 
        std::vector<fp_t> m1w0,
        std::vector<fp_t> m0w1,
        int record,
        int64_t tr) 
{
    size_t k = 0;
    for (uint i : operands) {
        if (tr >= 0) {
            uint64_t rand = rng();
            stim::bit_ref x = x_table[i][tr],
                            z = z_table[i][tr],
                            l = leak_table[i][tr],
                            lock = lock_table[i][tr];
            bool rx = rand & 0x1,
                 rz = rand & 0x2;
            __measure(x, z, l, lock, rx, rz);
        } else {
            const size_t width = z_table[i].num_bits_padded();
            stim::simd_bits<SIMD_WIDTH> rand_x = stim::simd_bits<SIMD_WIDTH>::random(width, rng),
                                        rand_z = stim::simd_bits<SIMD_WIDTH>::random(width, rng);
            std::array<stim::simd_bits_range_ref<SIMD_WIDTH>, 6> args{{
                x_table[i],
                z_table[i],
                leak_table[i],
                lock_table[i],
                rand_x,
                rand_z
            }};
            for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>, 6> word_array) {
                auto &x = word_array[0],
                     &z = word_array[1],
                     &l = word_array[2],
                     &lock = word_array[3],
                     &rx = word_array[4],
                     &rz = word_array[5];
                __measure(x, z, l, lock, rx, rz);
            });
        }

        if (record >= 0) {
            if (tr >= 0) {
                // Do this manually.
                if (!lock_table[i][tr]) {
                    record_table[record][tr] = x_table[i][tr];
                    fp_t e = get_probability_sample_from_rng();
                    if (record_table[record][tr] == 0) {
                        if (e < m1w0[k]) record_table[record][tr] = 1;
                    } else {
                        if (e < m0w1[k]) record_table[record][tr] = 0;
                    }
                }
            } else {
                record_table[record].for_each_word(x_table[i], lock_table[i],
                        [] (auto& rec, auto& x, auto& lock)
                        {
                            rec = andnot(lock, x) | (rec & lock);
                        });
                error_channel_m(record, m1w0[k], m0w1[k], lock_table[i]);
            }
            record++;
            k++;
        }
    }
}

void
FrameSimulator::R(std::vector<uint> operands, int64_t tr) {
    for (uint i : operands) {
        if (tr >= 0) {
            uint64_t rand = rng();
            stim::bit_ref x = x_table[i][tr],
                            z = z_table[i][tr],
                            l = leak_table[i][tr],
                            lock = lock_table[i][tr];
            bool r = rand & 0x1;
            __reset(x, z, l, lock, r);
        } else {
            x_table[i] &= lock_table[i];
            leak_table[i] &= lock_table[i];

            stim::simd_bits rand = stim::simd_bits<SIMD_WIDTH>::random(z_table[i].num_bits_padded(), rng);
            z_table[i].for_each_word(rand, lock_table[i],
            [&] (auto& z, auto& r, auto& lock)
            {
                z = andnot(lock, r) | (z & lock);
            });
        }
    }
}

StateSimulator::label_1q_t
FrameSimulator::eDP1(uint q, uint64_t t) {
    auto p = rng() & 3;
    x_table[q][t] ^= (bool)(p & 1);
    z_table[q][t] ^= (bool)(p & 2);
    return IXYZ[p];
}

StateSimulator::label_1q_t
FrameSimulator::eX(uint q, uint64_t t) {
    x_table[q][t] ^= 1;
    return L_X;
}

StateSimulator::label_1q_t
FrameSimulator::eY(uint q, uint64_t t) {
    x_table[q][t] ^= 1;
    z_table[q][t] ^= 1;
    return L_Y;
}

StateSimulator::label_1q_t
FrameSimulator::eZ(uint q, uint64_t t) {
    z_table[q][t] ^= 1;
    return L_Z;
}

StateSimulator::label_1q_t
FrameSimulator::eL(uint q, uint64_t t) {
    leak_table[q][t] ^= 1;
    return L_L;
}

StateSimulator::label_2q_t 
FrameSimulator::eDP2(uint q1, uint q2, uint64_t t) {
    auto p = rng() & 15;
    x_table[q1][t] ^= (bool)(p & 1);
    z_table[q1][t] ^= (bool)(p & 2);
    x_table[q2][t] ^= (bool)(p & 4);
    z_table[q2][t] ^= (bool)(p & 8);
    return std::make_pair(IXYZ[p & 3], IXYZ[(p >> 2) & 3]);
}

StateSimulator::label_2q_t 
FrameSimulator::eLI(uint q1, uint q2, uint64_t t) {
    auto p = rng() % 3;
    bool c1 = (p == 0) || (p == 2);
    bool c2 = (p == 1) || (p == 2);
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
    return std::make_pair(IL[c1], IL[c2]);
}

#define QONTRA_USE_MESSY_LEAKAGE_TRANSPORT

StateSimulator::label_2q_t 
FrameSimulator::eLT(uint q1, uint q2, uint64_t t) {
#ifdef QONTRA_USE_MESSY_LEAKAGE_TRANSPORT
    bool c1 = leak_table[q2][t] & ~leak_table[q1][t];
    bool c2 = leak_table[q1][t] & ~leak_table[q2][t];
#else
    bool c1 = leak_table[q1][t] ^ leak_table[q2][t];
    bool c2 = c1;
#endif
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
    return std::make_pair(IL[c1], IL[c2]);
}

void
FrameSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(x_table);
    z_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(z_table);
    leak_table_cpy = stim::simd_bit_table<SIMD_WIDTH>(leak_table);
}

void
FrameSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    StateSimulator::rollback_where(pred);
    for (uint i = 0; i < n_qubits; i++) {
        copy_where(x_table_cpy[i], x_table[i], pred);
        copy_where(z_table_cpy[i], z_table[i], pred);
        copy_where(leak_table_cpy[i], leak_table[i], pred);
    }
}

}   // qontra
