/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

namespace qontra {

inline DetailedStimCircuit&
DetailedStimCircuit::operator=(const stim::Circuit& other) {
    stim::Circuit::operator=(other);
    return *this;
}

inline DetailedStimCircuit&
DetailedStimCircuit::operator=(const DetailedStimCircuit& other) {
    stim::Circuit::operator=(other);
    detection_event_to_color = other.detection_event_to_color;
    flag_detection_events = other.flag_detection_events;
    return *this;
}

inline void
DetailedStimCircuit::apply_errors(
        std::string name,
        std::vector<uint64_t> operands,
        std::vector<fp_t> errors,
        bool is_2q) 
{
    for (size_t i = 0; i < errors.size(); i++) {
        if (is_2q) {
            safe_append_ua(name, 
                    {static_cast<uint32_t>(operands[2*i]), static_cast<uint32_t>(operands[2*i+1])},
                    errors[i]);
        } else {
            safe_append_ua(name, {static_cast<uint32_t>(operands[i])}, errors[i]);
        }
    }
}

template <size_t W, size_t N, class FUNC> inline void
for_each_word(std::array<stim::simd_bits_range_ref<W>, N> ranges, FUNC body) {
    std::array<stim::bitword<W>, N> values;

    const size_t n = ranges[0].num_simd_words;
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < N; i++) values[i] = *(ranges[i].ptr_simd + j);
        body(values);
    }
}

inline bool
andnot(bool x, bool y) {
    return !x & y;
}

template <size_t W> inline stim::bitword<W>
andnot(stim::bitword<W> x, stim::bitword<W> y) {
    return x.andnot(y);
}

}   // qontra

namespace stim {

template <size_t W> inline bitword<W> operator~(bitword<W> w) {
    return w ^ bitword<W>::tile64(UINT64_MAX);
}

}   // stim
