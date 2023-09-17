/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "sim/frame_sim.h"

namespace qontra {

#define L_I     StateSimulator::label_1q_t::I
#define L_X     StateSimulator::label_1q_t::X
#define L_Y     StateSimulator::label_1q_t::Y
#define L_Z     StateSimulator::label_1q_t::Z
#define L_L     StateSimulator::label_1q_t::L

static const StateSimulator::label_1q_t IXYZ[4] = { L_I, L_X, L_Y, L_Z };
static const StateSimulator::label_1q_t IL[2] = { L_I, L_L };

void
FrameSimulator::H(std::vector<uint> operands) {
    for (uint i : operands) {
        x_table[i].for_each_word(
                z_table[i],
                lock_table[i],
                [&] (auto& x, auto& z, auto& lock)
                {
                    auto tmp = z;
                    z = (x & ~lock) | (z & lock);
                    x = (tmp & ~lock) | (x & lock);
                });
    }
}

void
FrameSimulator::X(std::vector<uint> operands) {
    for (uint i : operands) {
        // Flip X wherever lock = 0
        lock_table[i].invert_bits();
        x_table[i] ^= lock_table[i];
        lock_table[i].invert_bits();
    }
}

void
FrameSimulator::Z(std::vector<uint> operands) {
    for (uint i : operands) {
        // Flip Z wherever lock = 0
        lock_table[i].invert_bits();
        z_table[i] ^= lock_table[i];
        lock_table[i].invert_bits();
    }
}

void
FrameSimulator::S(std::vector<uint> operands) {}    // Do not implement.

void
FrameSimulator::CX(std::vector<uint> operands) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
    
        stim::simd_bits lock_both(lock_table[j1]);
        lock_both |= lock_table[j2];

        x_table[j1].for_each_word(
            z_table[j1],
            leak_table[j1],
            x_table[j2],
            z_table[j2],
            leak_table[j2],
            lock_both,
        [&] (auto& x1, auto& z1, auto& l1,
            auto& x2, auto& z2, auto& l2,
            auto& lock)
        {
                stim::simd_word rx1(rng(), rng());
                stim::simd_word rz1(rng(), rng());
                stim::simd_word rx2(rng(), rng());
                stim::simd_word rz2(rng(), rng());
                x1 ^= l2 & rx1 & ~lock;
                z1 ^= ((z2 & ~l2) | (l2 & rz1)) & ~lock;
                x2 ^= ((x1 & ~l1) | (l1 & rx2)) & ~lock;
                z2 ^= l1 & rz2 & ~lock;
        });
    }
}

void
FrameSimulator::M(
        std::vector<uint> operands, 
        std::vector<fp_t> m1w0,
        std::vector<fp_t> m0w1,
        int record) 
{
    uint opk = 0;
    for (uint i : operands) {
        x_table[i].for_each_word(leak_table[i], lock_table[i],
        [&] (auto& x, auto& l, auto& lock)
        {
            stim::simd_word r(rng(), rng());
            x = (((x & ~l) | (r & l)) & ~lock) | (x & lock);
        });

        stim::simd_bits rand = 
            stim::simd_bits::random(z_table[i].num_bits_padded(), rng);
        z_table[i].for_each_word(rand, lock_table[i],
        [&] (auto& z, auto& r, auto& lock)
        {
            z = (r & ~lock) | (z & lock);
        });

        if (record >= 0) {
            record_table[record].for_each_word(x_table[i], lock_table[i],
                    [&] (auto& rec, auto& x, auto& lock)
                    {
                        rec = (x & ~lock) | (rec & lock);
                    });
            error_channel_m(record, m1w0[opk], m0w1[opk], lock_table[i]);
            record++;
            opk++;
        }
    }
}

void
FrameSimulator::R(std::vector<uint> operands) {
    for (uint i : operands) {
        x_table[i] &= lock_table[i];
        leak_table[i] &= lock_table[i];

        stim::simd_bits rand = 
            stim::simd_bits::random(z_table[i].num_bits_padded(), rng);
        z_table[i].for_each_word(rand, lock_table[i],
        [&] (auto& z, auto& r, auto& lock)
        {
            z = (r & ~lock) | (z & lock);
        });
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

StateSimulator::label_2q_t 
FrameSimulator::eLT(uint q1, uint q2, uint64_t t) {
    bool c1 = leak_table[q2][t] & ~leak_table[q1][t];
    bool c2 = leak_table[q1][t] & ~leak_table[q2][t];
    leak_table[q1][t] ^= c1;
    leak_table[q2][t] ^= c2;
    return std::make_pair(IL[c1], IL[c2]);
}

void
FrameSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = stim::simd_bit_table(x_table);
    z_table_cpy = stim::simd_bit_table(z_table);
    leak_table_cpy = stim::simd_bit_table(leak_table);
}

void
FrameSimulator::rollback_where(stim::simd_bits_range_ref pred) {
    StateSimulator::rollback_where(pred);
    for (uint i = 0; i < n_qubits; i++) {
        copy_where(x_table_cpy[i], x_table[i], pred);
        copy_where(z_table_cpy[i], z_table[i], pred);
        copy_where(leak_table_cpy[i], leak_table[i], pred);
    }
}

}   // qontra
