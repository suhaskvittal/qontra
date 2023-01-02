/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#include "fleece.h"

namespace qrc {

#define FLEECE_DEBUG

static const uint64_t MAX_SHOTS = 100000;
static std::mt19937_64 LOCAL_RNG(0);

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
    meas_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    leak_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    lattice_graph(),
    path_table(),
    curr_min_detector(0),
    curr_max_detector(0),
    rng(rng)
{
    lattice_graph = fleece::to_lattice_graph(circuit);
    path_table = fleece::compute_path_table(lattice_graph);

    if (double_stabilizer) {
        curr_max_detector = detectors_per_round >> 1;
    } else {
        curr_max_detector = detectors_per_round;
    }

#ifdef FLEECE_DEBUG
    // Print out lattice graph.
    for (auto v : lattice_graph.vertices()) {
        std::cout << v->qubit << " | is_data=" 
            << v->is_data << ", detectors={";
        for (auto d : v->detectors) {
            std::cout << " " << d;
        }
        std::cout << " }, mtimes={";
        for (auto t : v->measurement_times) {
            std::cout << " " << t;
        }
        std::cout << " }\n";
        std::cout << "\tadj={";
        for (auto w : lattice_graph.adjacency_list(v)) {
            std::cout << " " << w->qubit;
        }
        std::cout << " }\n";
    }
#endif
}

Fleece::SyndromeOutput
Fleece::generate_syndromes(uint64_t shots, uint last_leakage_round, uint64_t seed) {
    uint32_t row_size = circuit.count_detectors() + circuit.count_observables();
    meas_results = stim::simd_bit_table(row_size, shots);
    leak_results = stim::simd_bit_table(row_size, shots);

    auto det_obs = stim::DetectorsAndObservables(circuit);

    rng.seed(seed);
    sim.reset_all();
    curr_min_detector = 0;
    if (double_stabilizer) {
        curr_max_detector = detectors_per_round >> 1;
    } else {
        curr_max_detector = detectors_per_round;
    }

    uint round = 0;
    while (!sim.cycle_level_simulation(circuit)) {
        meas_results.clear();
        leak_results.clear();
        stim::read_from_sim(sim, det_obs, false, false, true, 
                meas_results, leak_results, curr_max_detector);
        if (!fake_run) {
            correct_leak(shots);
        }
        // Update curr min and curr max detectors.
        curr_min_detector += curr_max_detector; 
        curr_max_detector += detectors_per_round;
        if (curr_max_detector > circuit.count_detectors()) {
            curr_max_detector = circuit.count_detectors();
        }

        round++;
        if (round == last_leakage_round) {
            sim.toggle_leakage();
        }
    }

    meas_results.clear();
    leak_results.clear();
    stim::read_from_sim(sim, det_obs, false, true, true, meas_results, leak_results);
    SyndromeOutput out = std::make_pair(meas_results, leak_results);
    return out;
}

void
Fleece::correct_leak(uint64_t shots) {
    meas_results = meas_results.transposed();
    leak_results = leak_results.transposed();
    for (uint64_t s = 0; s < shots; s++) {
        // We only care about shots where there was a leakage.
        if (!leak_results[s].not_zero()) {
            continue;
        }
        
        for (uint32_t i = curr_min_detector; i < curr_max_detector; i++) {
            if (leak_results[s][i]) {
                auto vleak = lattice_graph.get_vertex_by_detector(i);
                if (vleak->is_data) {
                    correct_data_qubit(s, i, vleak);
                } else {
                    correct_parity_qubit(s, i, vleak);
                }
            }
        }
    } 
    leak_results = leak_results.transposed();
    meas_results = meas_results.transposed();
}

void
Fleece::correct_data_qubit(uint64_t shot_number, uint detector,
        fleece::LatticeGraph::Vertex * vleak) 
{
    // The leakage (in the worst case) could have only started after the last
    // time this data qubit was used as a proxy parity qubit (if at all).
    uint start = base_detector(detector);
    for (uint d : vleak->detectors) {
        if (detector == d) {
            break;
        }
        start = d;
    }
#ifdef FLEECE_DEBUG
    std::cout << "Correcting leak in shot=" << shot_number 
        << ", detector=" << detector << " from qubit " << vleak->qubit 
        << " (start=" << start << ")\n";
#endif
    // Iteratively, clear any flips in potentially affected rounds and
    // adjacent parity qubits to vleak.
    while (start <= detector) {
        const uint level = get_level(start);
        // Clear errors on prior rounds.
        meas_results[shot_number][start] = 0;
        const uint32_t mtime = get_measurement_time(start);
        if (level == 0) {
            sim.m_record.storage[mtime][shot_number] = 0;
        } else {
            const uint32_t pmtime = get_measurement_time(prev_detector(start));
            sim.m_record.storage[mtime][shot_number] = 
                sim.m_record.storage[pmtime][shot_number];
        }
        // Clear errors in adjacent parity qubits.
        for (auto par : lattice_graph.adjacency_list(vleak)) {
            if (par->detectors[0] == -1) {
                continue;
            }
            const uint d = jump_to_level(par->detectors[0], level);
            if (d == BOUNDARY_INDEX) {
                continue;
            }
            const uint32_t t = get_measurement_time(d);
#ifdef FLEECE_DEBUG
            std::cout << "\tUpdating measurement record for time=" << t
                << " for neighbor " << par->qubit << "(d=" << d << "):"
                << " was " << sim.m_record.storage[t][shot_number]
                << " (detector = " << meas_results[shot_number][d] << ")\n";
#endif
            meas_results[shot_number][d] = 0;
            if (level == 0) {
                sim.m_record.storage[t][shot_number] = 0;
            } else {
                const uint32_t pt = get_measurement_time(prev_detector(d));
                std::cout << "\t\tusing measurement at time " << pt << "\n";
                sim.m_record.storage[t][shot_number] = 
                    sim.m_record.storage[pt][shot_number];
            }
#ifdef FLEECE_DEBUG
            std::cout << "\t\tchanged to " 
                << sim.m_record.storage[t][shot_number] << "\n";
#endif
        } 

        start = next_detector(start);
    }
    // Clear the leakage.
    sim.leak_record.storage[get_measurement_time(detector)][shot_number] = 0;
}

void
Fleece::correct_parity_qubit(uint64_t shot_number, uint detector,
        fleece::LatticeGraph::Vertex * vleak) 
{

}

uint
Fleece::get_level(uint detector) {
    return (detector + (detectors_per_round>>1)) / detectors_per_round;
}

uint
Fleece::jump_to_level(uint base, uint level) {
    const uint32_t offset = base >= (detectors_per_round>>1);
    if (level == 0 && offset > 0) {
        return BOUNDARY_INDEX;
    } else if (level-offset == 0) {
        return base;
    } else {
        return base + (level-offset)*detectors_per_round;
    }
}

uint
Fleece::base_detector(uint detector) {
    return detector % detectors_per_round;
}

uint
Fleece::prev_detector(uint detector) {
    return detector - detectors_per_round;
}

uint
Fleece::next_detector(uint detector) {
    return detector + detectors_per_round;
}

uint
Fleece::get_measurement_time(uint detector) {
    const uint level = get_level(detector);
    const uint base = base_detector(detector);
    auto v = lattice_graph.get_vertex_by_detector(base);
    return v->measurement_times[0] + detectors_per_round * level;
}

stim::simd_bit_table
double_stabilizer_to_single_stabilizer(
    stim::simd_bit_table orig,
    uint code_dist,
    uint num_detectors,
    uint num_observables,
    uint num_shots,
    bool is_memory_z) 
{
    stim::simd_bit_table mod(num_detectors + num_observables, num_shots);

    const uint detectors_per_round = code_dist*code_dist - 1;
    const uint dprdiv2 = detectors_per_round >> 1;
    uint32_t orig_ptr = 0, mod_ptr = 0;
    while (mod_ptr < num_detectors) {
        if (orig_ptr % detectors_per_round < dprdiv2
            || mod_ptr >= num_detectors - dprdiv2) 
        {
            mod[mod_ptr++] |= orig[orig_ptr++];
        } else {
            orig_ptr++;
    void tower_correct(uint64_t shots);
    void clean_parity_qubit(uint32_t, uint lower_level, uint upper_level, uint64_t shot);
        }
    }
    for (uint i = 0; i < num_observables; i++) {
        mod[mod_ptr++] |= orig[orig_ptr++];
    }
    return mod;
}

} // qrc
