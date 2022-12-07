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
                char output_basis,
                uint32_t tower_cutoff)
    :circuit(circuit),
    detectors_per_round(detectors_per_round),
    double_stabilizer(double_stabilizer),
    reset_basis(reset_basis),
    output_basis(output_basis),
    tower_cutoff(tower_cutoff),
    sim(circuit.count_qubits(), MAX_SHOTS, SIZE_MAX, rng),
    meas_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    leak_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    lattice_graph(),
    curr_min_detector(0),
    curr_max_detector(0),
    rng(rng)
{
    lattice_graph = fleece::to_lattice_graph(circuit);

    if (double_stabilizer) {
        curr_max_detector = detectors_per_round >> 1;
    } else {
        curr_max_detector = detectors_per_round;
    }
}

Fleece::SyndromeOutput
Fleece::generate_syndromes(uint64_t shots) {
    uint32_t row_size = circuit.count_detectors() + circuit.count_observables();
    meas_results = stim::simd_bit_table(row_size, shots);
    leak_results = stim::simd_bit_table(row_size, shots);

    auto det_obs = stim::DetectorsAndObservables(circuit);
    while (!sim.cycle_level_simulation(circuit)) {
        stim::read_from_sim(sim, det_obs, false, true, true, meas_results, leak_results);
        // Correct data qubit "tower-like" errors.
        tower_correct(shots);


        // Update curr min and curr max detectors.
        curr_min_detector = curr_max_detector;
        curr_max_detector += detectors_per_round;
        if (curr_max_detector > circuit.count_detectors()) {
            curr_max_detector = circuit.count_detectors();
        }
    }

    stim::read_from_sim(sim, det_obs, false, true, false, meas_results, leak_results);
    SyndromeOutput out = std::make_pair(meas_results, leak_results);
    return out;
}


void
Fleece::tower_correct(uint64_t shots) {
    // First detect which data qubits have leaked. We do this by identifying 
    // all "towers" syndrome bits and group each tower based on adjacent data
    // qubits. Groups of towers with >1 syndrome bits for a data qubit are 
    // cleared and the FTQC is directed to reset the corresponding data qubits.
    uint detectors = curr_max_detector - curr_min_detector;
    
    std::vector<std::vector<uint8_t>> table(
            shots, std::vector<uint8_t>(detectors, 0));
    for (uint i = curr_min_detector; i < curr_max_detector; i++) {
        uint d = base_detector(i);
        for (uint64_t s = 0; s < shots; s++) {
            table[s][d] += meas_results[i][s];
        }
    }

    for (uint64_t s = 0; s < shots; s++) {
        std::map<uint32_t, std::vector<uint32_t>> marking_table;
        std::vector<uint8_t> freq_row(table[s]);
        for (uint i = 0; i < detectors; i++) {
            if (freq_row[i] >= tower_cutoff) {
                auto v = lattice_graph.get_vertex_by_detector(i);
                auto adjlist = lattice_graph.get_adjacency_list(v);
                for (const auto& w : adjlist) {
                    if (!marking_table.count(w.qubit)) {
                        marking_table[w.qubit] = std::vector<uint32_t>();
                    }
                    marking_table[w.qubit].push_back(v.qubit);
                }
            }
        }

        // Collect all the parity qubits that have
        // been affected by a leakage.
        std::set<uint32_t> clear_parity_qubits;
        for (auto entry : marking_table) {
            if (entry.second.size() > 1) {
                for (uint32_t parity_qubit : entry.second) {
                    clear_parity_qubits.insert(parity_qubit);
                }
                // Also, reset the data qubit.
                uint32_t q = entry.first;
                sim.x_table[q][s] = 0;
                if (sim.guarantee_anticommutation_via_frame_randomization) {
                    sim.z_table[q][s] = rng() & 0x1;
                }
                sim.leakage_table[q][s] = 0;
            }
        }

        for (uint32_t q : clear_parity_qubits) {
            // Clear bits in the measurement and leakage records
            // and the syndromes.
            const auto& v = lattice_graph.get_vertex_by_qubit(q);             
            uint currd = v.base_detector;
            while (currd < curr_max_detector) {
                meas_results[s][currd] = 0;
                leak_results[s][currd] = 0;
                currd = next_detector(currd);
            }

            for (uint32_t meas_t : v.measurement_times) {
                sim.m_record.storage[meas_t][s] = 0; 
                sim.leak_record.storage[meas_t][s] = 0; 
            }
        }
    }
}

uint
Fleece::base_detector(uint detector) {
    /*
    if (double_stabilizer) {
        const uint dprdiv2 = detectors_per_round >> 1;
        if (detector < dprdiv2) {
            return detector;
        } else if (detector >= circuit.count_detectors() - dprdiv2) {
            return detector - (circuit.count_detectors() - dprdiv2);
        } else {
            return detector % detectors_per_round;
        }
    } else {
        return detector % detectors_per_round;
    }
    */
    return detector % detectors_per_round;
}

uint
Fleece::next_detector(uint detector) {
    return detector + detectors_per_round;
}



} // qrc
