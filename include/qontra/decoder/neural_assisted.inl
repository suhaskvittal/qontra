/*
 *  author: Suhas Vittal
 *  date:   19 Feburary 2024
 * */

namespace qontra {

inline std::vector<uint64_t>
NeuralAssistedDecoder::get_flags(std::vector<uint64_t>& detectors) {
    std::vector<uint64_t> flags;
    for (auto it = detectors.begin(); it != detectors.end(); ) {
        if (circuit.flag_detectors.count(*it)) {
            size_t k = 0;
            for (uint64_t f : circuit.flag_detectors) {
                if (f == *it) {
                    flags.push_back(k);
                }
                k++;
            }
            it = detectors.erase(it);
        } else {
            it++;
        }
    }
    return flags;
}

}   // qontra
