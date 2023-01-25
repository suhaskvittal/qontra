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
#ifdef FLEECE_DEBUG
            std::cout << "round = " << round << "\n";
#endif
            correct_leak(shots);
        }
        // Update curr min and curr max detectors.
        curr_min_detector = curr_max_detector; 
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

std::map<uint32_t, uint8_t>
Fleece::get_data_qubit_states(uint64_t shot_number) {
    std::map<uint32_t, uint8_t> state_table; 
    for (auto v : lattice_graph.vertices()) {
        if (v->is_data) {
            uint8_t x = sim.x_table[v->qubit][shot_number];
            uint8_t z = sim.z_table[v->qubit][shot_number];
            state_table[v->qubit] = (x << 1) | z;
        }
    }
    return state_table;
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
                if (await_detector_set.count(i)) {
                    uint32_t j = prev_detector(i);
                    auto vleak = lattice_graph.get_vertex_by_detector(j);
                    // We need the results of the next round to accurately
                    // kill leakage errors.
                    if (vleak->is_data) {
                        correct_data_qubit(s, i, vleak);
                    } else {
                        correct_parity_qubit(s, i, vleak);
                    }
                } else if (next_detector(i) < 
                        circuit.count_detectors() - (detectors_per_round >> 1))
                {
                    await_detector_set.insert(next_detector(i));
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
        if (prev_detector(detector) == d) {
            break;
        }
        start = d;
    }
#ifdef FLEECE_DEBUG
    std::cout << "[data] Correcting leak in shot=" << shot_number 
        << ", detector=" << detector << " from qubit " << vleak->qubit 
        << " (start=" << start << ")\n";
#endif
    // Iteratively, search for flips in detectors adjacent to the leaked qubit.
    // If we hit a critical point, reset measurements in those detectors.
    const uint level = get_level(start);
    std::vector<uint> adjacent_detectors;
#ifdef FLEECE_DEBUG
    std::cout << "\tadjacent detectors:";
#endif
    for (auto w : lattice_graph.adjacency_list(vleak)) {
        if (w->detectors[0] == -1) {
            continue;
        }
        uint base = base_detector(w->detectors[0]);
        adjacent_detectors.push_back(jump_to_level(base, level));
#ifdef FLEECE_DEBUG
        std::cout << " " << base;
#endif
    }
#ifdef FLEECE_DEBUG
    std::cout << "\n";
#endif
    uint curr_level = level;
    const uint target_level = get_level(detector);
    while (curr_level <= target_level) {
        // Clear measurements
        for (uint d : adjacent_detectors) {
            if (get_level(d) != curr_level) {
                continue;
            }
            meas_results[shot_number][d] = 0;
            const uint32_t mt = get_measurement_time(d);
            if (curr_level == 0) {
                sim.m_record.storage[mt][shot_number] = 0;
            } else {
                const uint32_t pmt = get_measurement_time(prev_detector(d));
                sim.m_record.storage[mt][shot_number] = 
                    sim.m_record.storage[pmt][shot_number];
            }
        }

        for (uint i = 0; i < adjacent_detectors.size(); i++) {
            if (get_level(adjacent_detectors[i]) != curr_level) {
                continue;
            }
            adjacent_detectors[i] = next_detector(adjacent_detectors[i]);
        }
        curr_level++;
    }
    // Clear the leakage.
    const uint32_t mt = get_measurement_time(detector);
    const uint32_t pmt = get_measurement_time(prev_detector(detector));
    sim.leak_record.storage[mt][shot_number] = 0;
    sim.leak_record.storage[pmt][shot_number] = 0;
}

void
Fleece::correct_parity_qubit(uint64_t shot_number, uint detector,
        fleece::LatticeGraph::Vertex * vleak) 
{
    // The leakage (in the worst case) could have only started after the last
    // time this parity qubit was measured.
    uint start = base_detector(detector);
    for (uint d : vleak->detectors) {
        if (prev_detector(detector) == d) {
            break;
        }
        start = d;
    }
#ifdef FLEECE_DEBUG
    std::cout << "[parity] Correcting leak in shot=" << shot_number 
        << ", detector=" << detector << " from qubit " << vleak->qubit 
        << " (start=" << start << ")\n";
#endif
    uint curr_level = get_level(start); 
    const uint max_level = get_level(detector);
    
    std::map<uint32_t, uint8_t> paradjust;

    while (curr_level <= max_level) {
        // Find other detectors which are active.
        std::vector<uint> other_active;
        uint ds = curr_level * detectors_per_round - (detectors_per_round>>1);
        uint dt;
        if (curr_level == 0) {
            dt = detectors_per_round >> 1;
        } else {
            dt = ds + detectors_per_round;
        }
        for (uint i = ds; i < dt; i++) {
            if (!paradjust.count(base_detector(i))) {
                paradjust[base_detector(i)] = 0;
            }
            const uint8_t offset = paradjust[base_detector(i)];
            if (meas_results[shot_number][i] ^ offset) {
                other_active.push_back(i);
            }
        }
#ifdef FLEECE_DEBUG
        std::cout << "\tclearing active detectors:\n";
#endif
        std::set<fleece::LatticeGraph::Vertex*> xflipped;
        std::set<fleece::LatticeGraph::Vertex*> zflipped;
        for (uint d : other_active) {
            auto vbase = lattice_graph.get_vertex_by_detector(base_detector(d));
            auto common = lattice_graph.get_common_neighbors(vleak, vbase);
            if (common.empty()) {
                continue;
            }
#ifdef FLEECE_DEBUG
            std::cout << "\t" << d << "\treset:";
#endif
            // If there are shared data qubits, flip them.
            // Unless this detector has the same base (first round) detector
            // as the leaked detector. Then, only clear the measurement.
            if (base_detector(d) != base_detector(detector)) {
                for (auto vdata : common) {
#ifdef FLEECE_DEBUG
                    std::cout << " " << vdata->qubit;
#endif
                    bool flippedx = false;
                    bool flippedz = false;
                    if (base_detector(d) < (detectors_per_round>>1)
                        && !xflipped.count(vdata)) 
                    {
#ifdef FLEECE_DEBUG
                        std::cout << "(X)";
#endif
                        sim.x_table[vdata->qubit][shot_number] ^= 1;
                        xflipped.insert(vdata);
                        flippedx = true;
                    } 
                    if (base_detector(d) >= (detectors_per_round>>1)
                        && !zflipped.count(vdata))
                    {
#ifdef FLEECE_DEBUG
                        std::cout << "(Z)";
#endif
                        sim.z_table[vdata->qubit][shot_number] ^= 1;
                        zflipped.insert(vdata);
                        flippedz = true;
                    }
                    // Update parity adjustments.
                    for (auto w : lattice_graph.adjacency_list(vdata)) {
                        if (w->detectors[0] == -1) {
                            continue;
                        }
                        if (!paradjust.count(w->detectors[0])) {
                            paradjust[w->detectors[0]] = 0;
                        }
                        uint8_t update;
                        if (w->detectors[0] < (detectors_per_round>>1)) {
                            update = flippedx;
                        } else {
                            update = flippedz;
                        }
                        paradjust[w->detectors[0]] ^= update;
                    }
                }
            }
            // Then, clear the measurements for the detector.
            const uint32_t mt = get_measurement_time(d);
            if (!paradjust.count(base_detector(d))) {
                paradjust[base_detector(d)] = 0;
            }
            const uint8_t offset = paradjust[base_detector(d)];
            if (curr_level == 0) {
                sim.m_record.storage[mt][shot_number] = 0;
            } else {
                const uint32_t pmt = get_measurement_time(prev_detector(d));
                sim.m_record.storage[mt][shot_number] = 
                    sim.m_record.storage[pmt][shot_number];
            }
#ifdef FLEECE_DEBUG
            std::cout << "\n";
#endif
        }
        curr_level++;
    }
    // Clear leakage
    const uint32_t mt = get_measurement_time(detector);
    const uint32_t pmt = get_measurement_time(prev_detector(detector));
    sim.leak_record.storage[mt][shot_number] = 0;
    sim.leak_record.storage[pmt][shot_number] = 0;
}

uint
Fleece::get_level(uint detector) {
    return (detector + (detectors_per_round>>1)) / detectors_per_round;
}

uint
Fleece::jump_to_level(uint base, uint level) {
    const uint32_t offset = base >= (detectors_per_round>>1);
    if (level == 0 && offset > 0) {
        return base;
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
