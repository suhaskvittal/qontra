/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "sim/frame_sim.h"

namespace qontra {

void
FrameSimulator::H(std::vector<uint> operands) {
    for (uint i : operands) x_table[i].swap_with(z_table[i]);
}

void
FrameSimulator::X(std::vector<uint> operands) {
    for (uint i : operands) x_table[i].invert_bits();
}

void
FrameSimulator::Z(std::vector<uint> operands) {
    for (uint i : operands) z_table[i].invert_bits();
}

void
FrameSimulator::S(std::vector<uint> operands) {}    // Do not implement.

void
FrameSimulator::CX(std::vector<uint> operands) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
        x_table[j1].for_each_word(
            z_table[j1],
            leak_table[j1],
            x_table[j2],
            z_table[j2],
            leak_table[j2],
        [&] (auto& x1, auto& z1, auto& l1, auto& x2, auto& z2, auto& l2)
        {
                stim::simd_word rx1(rng(), rng());
                stim::simd_word rz1(rng(), rng());
                stim::simd_word rx2(rng(), rng());
                stim::simd_word rz2(rng(), rng());
                x1 ^= l2 & rx1;
                z1 ^= (z2 & ~l2) | (l2 & rz1);
                x2 ^= (x1 & ~l1) | (l1 & rx2);
                z2 ^= l1 & rz2;
        });
    }
}

void
FrameSimulator::M(std::vector<uint> operands, bool record) {
    for (uint i : operands) {
        x_table[i].for_each_word(leak_table[i], 
        [&] (auto& x, auto& l)
        {
            stim::simd_word r(rng(), rng());
            x = (x & ~l) | (r & l);
        });
        z_table[i].randomize(z_table[i].num_bits_padded(), rng);
        if (record) {
            record_table[record_offset].clear();
            record_table[record_offset++] |= x_table[i];
        }
    }
}

void
FrameSimulator::R(std::vector<uint> operands) {
    for (uint i : operands) {
        x_table[i].clear();
        z_table[i].randomize(z_table[i].num_bits_padded(), rng);
        leak_table[i].clear();
    }
}

void
FrameSimulator::eDPO(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            // Implement as reset. We conservatively assume that
            // if the qubit was in |2>, it has only relaxed to |1>.
            x_table[j][t] = leak_table[j][t];
            z_table[j][t] = rng();
            leak_table[j][t] = 0;
        });
    }
}

void
FrameSimulator::eDPH(std::vector<uint> operands, std::vector<fp_t> rates) {
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
            auto tmp = x_table[j][t];
            x_table[j][t] = z_table[j][t];
            z_table[j][t] = tmp;
        });
    }
}

void
FrameSimulator::eDP1(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            auto p = rng() & 3;
            x_table[j][t] ^= (p & 1);
            z_table[j][t] ^= (p & 2);
        });
    }
}

void
FrameSimulator::eDP2(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) {
            auto p = rng() & 15;
            x_table[j1][t] ^= (p & 1);
            z_table[j1][t] ^= (p & 2);
            x_table[j2][t] ^= (p & 4);
            z_table[j2][t] ^= (p & 8);
        });
    }
}

void
FrameSimulator::eX(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            x_table[j][t] ^= 1;
        });
    }
}

void
FrameSimulator::eLI(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (uint i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i], j2 = operands[i+1];
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
FrameSimulator::eLT(std::vector<uint> operands, std::vector<fp_t> rates) {
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
FrameSimulator::snapshot() {
    StateSimulator::snapshot();
    x_table_cpy = stim::simd_bit_table(x_table);
    z_table_cpy = stim::simd_bit_table(z_table);
    leak_table_cpy = stim::simd_bit_table(leak_table);
}

void
FrameSimulator::rollback_at_trial(uint64_t t) {
    StateSimulator::rollback_at_trial(t);
    for (uint i = 0; i < n_qubits; i++) {
        x_table[i][t] = x_table_cpy[i][t];
        z_table[i][t] = z_table_cpy[i][t];
        leak_table[i][t] = leak_table_cpy[i][t];
    }
}

}   // qontra
