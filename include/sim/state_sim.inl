/* author: Suhas Vittal
 *  date:   25 December 2023
 * */

namespace qontra {

inline void
StateSimulator::set_seed(uint64_t x) {
    rng.seed(x);
}

inline void
StateSimulator::reset_sim() {
    record_table.clear();
    lock_table.clear();
}

template <StateSimulator::ErrorChannel1Q CH> inline void
StateSimulator::error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (size_t i = 0; i < operands.size(); i++) {
        uint j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            if (lock_table[j][t])   return;
            (this->*CH)(j, t);
        });
    }
}

template <StateSimulator::ErrorChannel2Q CH> inline void
StateSimulator::error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
    for (size_t i = 0; i < operands.size(); i += 2) {
        uint j1 = operands[i];
        uint j2 = operands[i+1];
        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) 
        {
            if (lock_table[j1][t] || lock_table[j2][t]) return;
            (this->*CH)(j1, j2, t);
        });
    }
}

void
StateSimulator::error_channel_m(
        uint64_t rec,
        fp_t m1w0,
        fp_t m0w1,
        stim::simd_bits_range_ref<SIMD_WIDTH> lock) 
{
    stim::RareErrorIterator::for_samples(m1w0, shots, rng,
            [&] (size_t t)
            {
                if (lock[t]) return;
                if (record_table[rec][t] == 0) record_table[rec][t] = 1;
            });
    stim::RareErrorIterator::for_samples(m0w1, shots, rng,
            [&] (size_t t)
            {
                if (lock[t]) return;
                if (record_table[rec][t] == 1) record_table[rec][t] = 0;
            });
}

inline void
StateSimulator::shift_record_by(uint64_t offset) {
    for (size_t i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        if (i < offset) record_table[i].clear();
        else            record_table[i].swap_with(record_table[i-offset]);
    }
}

inline void
StateSimulator::snapshot() {
    record_table_cpy = stim::simd_bit_table(record_table);
    lock_table_cpy = stim::simd_bit_table(lock_table);
}

inline void
StateSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    for (size_t i = 0; i < statesim::G_RECORD_SPACE_SIZE; i++) {
        copy_where(record_table_cpy[i], record_table[i], pred);
    }
    for (size_t i = 0; i < n_qubits; i++) {
        copy_where(lock_table_cpy[i], lock_table[i], pred);
    }
}

}   // qontra
