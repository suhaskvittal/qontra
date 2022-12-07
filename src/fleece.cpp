/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#include "fleece.h"

namespace qrc {

Fleece::Fleece(const stim::Circuit& circuit,
                std::mt19937_64& rng,
                uint detectors_per_round,
                bool double_stabilizer,
                char reset_basis,
                char output_basis)
    :circuit(circuit),
    detectors_per_round(detectors_per_round),
    double_stabilizer(double_stabilizer),
    reset_basis(reset_basis),
    output_basis(output_basis),
    sim(circuit.count_qubits(), MAX_SHOTS, SIZE_MAX, rng),
    current_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    decoding_graph(),
    path_table(),
    stabilizer_to_data(),
    curr_min_detector(0)
{
    decoding_graph = to_decoding_graph(circuit);
    path_table = compute_path_table(decoding_graph);

    // Traverse the circuit to populate stabilizer_to_data.
    for (auto op : circuit.operations) {
        std::string name(op->gate.name);
        if (name != "CX" && name != "ZCX") {
            continue;
        }
        auto& targets = op.target_data.targets;
        for (uint i = 0; i < targets.size(); i++) {
             
        }
    }
}


void
Fleece::tower_correct(uint64_t shots) {
    const curr_max_detector = curr_min_detector + detectors_per_round;

    // First detect which data qubits have leaked. We do this by identifying 
    // all "towers" syndrome bits and group each tower based on adjacent data
    // qubits. Groups of towers with >1 syndrome bits for a data qubit are 
    // cleared and the FTQC is directed to reset the corresponding data qubits.
    std::vector<std::vector<uint8_t>> table(
            shots, std::vector<uint8_t>(detectors_per_round, 0));
    for (uint i = 0; i < curr_max_detector; i++) {
        uint d = base_detector(i);
        for (uint64_t s = 0; s < shots; s++) {
            table[s][d] += current_results[i][s];
        }
    }

    
}

uint
Fleece::base_detector(uint detector) {
    if (double_stabilizer) {
        const dprdiv2 = detectors_per_round >> 1;
        if (detector < dprdiv2) {
            return detector;
        } else if (detector >= circuit.count_detectors() - dprdiv2) {
            return detector - (circuit.count_detectors() - dprdiv2);
        } else {
            return (detector + dprdiv2) % detectors_per_round;
        }
    } else {
        return detector % detectors_per_round;
    }
}


} // qrc
