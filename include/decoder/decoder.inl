/*
 *  author: Suhas Vittal
 *  date    25 December 2023
 * */

inline DetailedStimCircuit
Decoder::get_circuit() {
    return circuit;
}

template <class T> inline std::vector<uint>
Decoder::get_nonzero_detectors_(T syndrome) {
    return get_nonzero_detectors(syndrome, circuit.count_detectors());
}
