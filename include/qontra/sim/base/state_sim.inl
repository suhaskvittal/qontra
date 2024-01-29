/* 
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

namespace qontra {

template <class T, class U> inline void
copy_where(T from, T to, U pred) {
    from.for_each_word(to, pred, 
            [&] (auto& f, auto& t, auto& p)
            {
                t = (p & f) | p.andnot(t);
            });
}

// Specialization to T = stim::simd_bit_table<SIMD_WIDTH>&
template <class U> inline void
copy_where(stim::simd_bit_table<SIMD_WIDTH>& from, stim::simd_bit_table<SIMD_WIDTH>& to, U pred) {
    for (size_t i = 0; i < from.num_major_bits_padded(); i++) {
        copy_where(from[i], to[i], pred);
    }
}

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
StateSimulator::error_channel(std::vector<uint64_t> operands, std::vector<fp_t> rates) {
    for (size_t i = 0; i < operands.size(); i++) {
        uint64_t j = operands[i];
        stim::RareErrorIterator::for_samples(rates[i], shots, rng,
        [&] (size_t t) {
            if (lock_table[j][t])   return;
            (this->*CH)(j, t);
        });
    }
}

template <StateSimulator::ErrorChannel2Q CH> inline void
StateSimulator::error_channel(std::vector<uint64_t> operands, std::vector<fp_t> rates) {
    for (size_t i = 0; i < operands.size(); i += 2) {
        uint64_t j1 = operands[i];
        uint64_t j2 = operands[i+1];
        stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
        [&] (size_t t) 
        {
            if (lock_table[j1][t] || lock_table[j2][t]) return;
            (this->*CH)(j1, j2, t);
        });
    }
}

inline void
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
StateSimulator::shift_record_by(int64_t offset) {
    left_shift(record_table, offset);
}

inline void
StateSimulator::snapshot() {
    record_table_cpy = stim::simd_bit_table(record_table);
    lock_table_cpy = stim::simd_bit_table(lock_table);
}

inline void
StateSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    copy_where(record_table_cpy, record_table, pred);
    copy_where(lock_table_cpy, lock_table, pred);
}

}   // qontra
