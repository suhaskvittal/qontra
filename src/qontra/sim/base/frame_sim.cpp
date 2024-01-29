/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "qontra/sim/base/frame_sim.h"
#include "qontra/sim/base/gate_impl/frame_sim.h"

namespace qontra {

FrameSimulator::FrameSimulator(uint64_t n, uint64_t max_shots)
    :StateSimulator(n, max_shots),
    x_table(n, max_shots),
    z_table(n, max_shots),
    leak_table(n, max_shots),
    x_table_cpy(n, max_shots),
    z_table_cpy(n, max_shots),
    leak_table_cpy(n, max_shots)
{}

FrameSimulator::FrameSimulator(const FrameSimulator& other)
    :StateSimulator(other),
    x_table(other.x_table),
    z_table(other.z_table),
    leak_table(other.leak_table),
    x_table_cpy(other.x_table_cpy),
    z_table_cpy(other.z_table_cpy),
    leak_table_cpy(other.leak_table_cpy)
{}

FrameSimulator::FrameSimulator(FrameSimulator&& other)
    :StateSimulator(other),
    x_table(std::move(other.x_table)),
    z_table(std::move(other.z_table)),
    leak_table(std::move(other.leak_table)),
    x_table_cpy(std::move(other.x_table_cpy)),
    z_table_cpy(std::move(other.z_table_cpy)),
    leak_table_cpy(std::move(other.leak_table_cpy))
{}

void
FrameSimulator::H(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i : operands) {
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
FrameSimulator::X(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i : operands) {
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
FrameSimulator::Z(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i : operands) {
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
FrameSimulator::S(std::vector<uint64_t> operands, int64_t tr) {} // Do not implement.

void
FrameSimulator::CX(std::vector<uint64_t> operands, int64_t tr) {
    for (size_t i = 0; i < operands.size(); i += 2) {
        uint64_t j1 = operands[i], j2 = operands[i+1];
    
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
            for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>*, 11>& word_array) 
            {
                auto& x1 = *word_array[0],
                    & z1 = *word_array[1],
                    & l1 = *word_array[2],
                    & x2 = *word_array[3],
                    & z2 = *word_array[4],
                    & l2 = *word_array[5],
                    & lock = *word_array[6],
                    & rx1 = *word_array[7],
                    & rz1 = *word_array[8],
                    & rx2 = *word_array[9],
                    & rz2 = *word_array[10];
                __cx_gate(x1, z1, l1, x2, z2, l2, lock, rx1, rz1, rx2, rz2);
            });
        }
    }
}

void
FrameSimulator::LEAKAGE_ISWAP(std::vector<uint64_t> operands, int64_t tr) {
    for (size_t i = 0; i < operands.size(); i += 2) {
        uint64_t j1 = operands[i], j2 = operands[i+1];
    
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
        std::vector<uint64_t> operands, 
        std::vector<fp_t> m1w0,
        std::vector<fp_t> m0w1,
        int record,
        int64_t tr) 
{
    size_t k = 0;
    for (uint64_t i : operands) {
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
            for_each_word(args, [] (std::array<stim::bitword<SIMD_WIDTH>*, 6> word_array) {
                auto &x = *word_array[0],
                     &z = *word_array[1],
                     &l = *word_array[2],
                     &lock = *word_array[3],
                     &rx = *word_array[4],
                     &rz = *word_array[5];
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
FrameSimulator::R(std::vector<uint64_t> operands, int64_t tr) {
    for (uint64_t i : operands) {
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

}   // qontra
