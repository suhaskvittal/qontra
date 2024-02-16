/*
 *  author: Suhas Vittal
 *  date    25 December 2023
 * */

#include <strings.h>

namespace qontra {

template <class T> std::vector<uint64_t>
get_nonzero_detectors_(T syndrome, uint64_t number_of_detectors) {
    std::vector<uint64_t> det;
    uint64_t w = 0;
    uint64_t last_bit = 0;
    while (det.size() < syndrome.popcnt()) {
        uint64_t i = ffsll(syndrome.u64[w] & ~((1L << last_bit)-1));
        if (i == 0) {   // No match found.
            last_bit = 0;
            w++;
            continue;
        }
        uint64_t d = (w << 6) | (i-1);
        if (d >= number_of_detectors) break;
        det.push_back(d);
        last_bit = i & 0x3f;
        w += (i >= 64);
    }
    return det;
}

inline DetailedStimCircuit
Decoder::get_circuit() {
    return circuit;
}

template <class T> inline std::vector<uint64_t>
Decoder::get_nonzero_detectors(T syndrome) {
    return get_nonzero_detectors_(syndrome, circuit.count_detectors());
}

}   // qontra
