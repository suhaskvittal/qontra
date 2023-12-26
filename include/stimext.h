/*
 *  author: Suhas Vittal
 *  date:   27 November 2023
 * */

#ifndef QONTRA_STIM_EXTENSIONS
#define QONTRA_STIM_EXTENSIONS

#include <stim.h>

#include <map>
#include <set>

namespace qontra {

// We provide an extension of Stim's Circuit class. The DetailedStimCircuit
// class can be used in any function accepting stim::Circuits, but has
// extra metadata that can be used by decoders and such.
//
// Stim's circuits cannot hold much information beyond error rates
// and operations, so we contain extra information in here.
struct DetailedStimCircuit : public stim::Circuit {
    DetailedStimCircuit() :stim::Circuit() {}
    DetailedStimCircuit(const stim::Circuit& other)
        :stim::Circuit(other)
    {}
    DetailedStimCircuit(const DetailedStimCircuit& other)
        :stim::Circuit(other),
        detection_event_to_color(other.detection_event_to_color),
        flag_detection_events(other.flag_detection_events),
        flag_edge_table(other.flag_edge_table)
    {}

    typedef std::tuple<uint, uint, uint>    flag_edge_t;

    std::map<uint, int> detection_event_to_color;
    std::set<uint>      flag_detection_events;

    std::map<uint, flag_edge_t> flag_edge_table;
};

}   // qontra

// Extension of stim::simd_bits_range_ref for_each_word function.
template <size_t W, size_t N, class FUNC> inline void
for_each_word(std::array<stim::simd_bits_range_ref<W>, N> ranges, FUNC body) {
    std::array<stim::bitword<W>, N> values;

    const size_t n = ranges[0].num_simd_words;
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < N; i++) values[i] = *(ranges[i].ptr_simd + j);
        body(values);
    }
}

// Template and specialization of ANDNOT for general use.
inline bool
andnot(bool x, bool y) { return ~x & y; }

template <size_t W> inline stim::bitword<W>
andnot(stim::bitword<W> x, stim::bitword<W> y) { return x.andnot(y); }

namespace stim {
    
template <size_t W> inline bitword<W> operator~(bitword<W> w) {
    return w ^ bitword<W>::tile64(UINT64_MAX);
}

}   // stim

#endif  // QONTRA_STIM_EXTENSIONS
