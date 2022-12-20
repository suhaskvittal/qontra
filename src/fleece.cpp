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
                char output_basis,
                uint32_t tower_cutoff,
                uint32_t sliding_window_size)
    :circuit(circuit),
    detectors_per_round(detectors_per_round),
    double_stabilizer(double_stabilizer),
    reset_basis(reset_basis),
    output_basis(output_basis),
    tower_cutoff(tower_cutoff),
    sliding_window_size(sliding_window_size),
    sim(circuit.count_qubits(), MAX_SHOTS, SIZE_MAX, rng),
    meas_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    leak_results(
            circuit.count_detectors() + circuit.count_observables(), MAX_SHOTS),
    lattice_graph(),
    path_table(),
    curr_min_detector(0),
    curr_max_detector(0),
    fake_run(false),
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
        if (v->base_detector >= 0) continue;
        std::cout << v->qubit << " | is_data=" 
            << v->is_data << ", base=" 
            << v->base_detector << ", mtimes={";
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
    // Print out path table
    for (auto v1 : lattice_graph.vertices()) {
        for (auto v2 : lattice_graph.vertices()) {
            std::cout << path_table[std::make_pair(v1, v2)].distance << "("
                << path_table[std::make_pair(v1, v2)].path.size() << ")"
                << " ";
        }
        std::cout << "\n";
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
#ifdef FLEECE_DEBUG
        std::cout << "(min, max) = (" << curr_min_detector << ", " << curr_max_detector << ")\n";
#endif
        stim::read_from_sim(sim, det_obs, false, true, true, meas_results, leak_results);
        // Correct data qubit "tower-like" errors.
        if (!fake_run) {
            tower_correct(shots);
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

bool
Fleece::toggle_fake_run() {
    fake_run = !fake_run;
    return fake_run;
}

#define TOWER_SIZE 30

void
Fleece::tower_correct(uint64_t shots) {
    // To make this fast, we need to optimize for the cache.
    // Thus, we transpose the meas_results and leak_results to improve cache hits.
    auto meas_results_T = meas_results.transposed();
    for (uint64_t s = 0; s < shots; s++) {
        // We first identify which data qubits have leaked by
        // checking the adjacent stabilizers of each data qubit.
        // If there is a tower larger than tower_cutoff, then
        // we declare a leak.
        if (!meas_results_T[s].not_zero()) {
            continue;
        }
        std::map<uint32_t, uint8_t> activity_table;
        std::vector<int32_t> leak_detections;
        for (uint i = 0; i < curr_max_detector; i++) {
            if (!meas_results_T[s][i]) {
                continue;
            }
            uint d = base_detector(i);
            if (!activity_table.count(d)) {
                activity_table[d] = 0;
            }
            if ((++activity_table[d]) == tower_cutoff) {
                leak_detections.push_back(d);
            }
        }
        if (leak_detections.size() & 0x1) {
            leak_detections.push_back(LATTICE_BOUNDARY_DETECTOR);
        }
        if (leak_detections.size() == 0) {
            continue;
        } else {
            std::cout << "leak detected with " 
                    << leak_detections.size() << " measurements in shot "
                    << s << ":";
            for (int32_t d : leak_detections) {
                if (d == LATTICE_BOUNDARY_DETECTOR) {
                    std::cout << " B";
                } else {
                    std::cout << " " << d;
                }
            }
            std::cout << "\n";
        }
        // Accumulate leakage detections and perform MWPM.
        const uint32_t n_mwpm_vertices = leak_detections.size();
        const uint32_t n_mwpm_edges = (n_mwpm_vertices * (n_mwpm_vertices-1)) >> 1;
        PerfectMatching pm(n_mwpm_vertices, n_mwpm_edges);

        for (uint i = 0; i < n_mwpm_vertices; i++) {
            int32_t di = leak_detections[i];
            auto vi = lattice_graph.get_vertex_by_detector(di);
            for (uint j = i+1; j < n_mwpm_vertices; j++) {
                int32_t dj = leak_detections[i];
                auto vj = lattice_graph.get_vertex_by_detector(dj);

                auto vi_vj = std::make_pair(vi, vj);
                pm.AddEdge(i, j, path_table[vi_vj].distance);
            }
        }
        pm.options.verbose = false;
        pm.Solve();

        // For each matched detector,
        // (1) clear all syndrome bits for those detectors for the suspected
        //  leaked rounds.
        // (2) reset the data qubits between the matched syndrome qubits and clear
        //  all parity qubit measurements in between the matched qubits.
        std::set<int32_t> visited;
        std::set<int32_t> affected_parity_qubits;
        for (uint i = 0; i < n_mwpm_vertices; i++) {
            int32_t di = leak_detections[i];
            if (visited.count(di)) {
                continue;
            }
            uint j = pm.GetMatch(i);
            int32_t dj = leak_detections[j];

            auto vi = lattice_graph.get_vertex_by_detector(di);
            auto vj = lattice_graph.get_vertex_by_detector(dj);
            auto vi_vj = std::make_pair(vi, vj);

            auto correcting_path = path_table[vi_vj].path;
            std::cout << "\tcorrecting(" << correcting_path.size() << "):";
            for (auto w : correcting_path) {
                int32_t q = w->qubit;
                std::cout << " " << q;
                if (w->is_data) {
                    std::cout << "(D)";
                    sim.x_table[q][s] = 0;
                    if (sim.guarantee_anticommutation_via_frame_randomization) {
                        sim.z_table[q][s] = LOCAL_RNG() & 0x1;
                    }
                    sim.leakage_table[q][s] = 0;
                } else {
                    if (q != LATTICE_BOUNDARY_DETECTOR) {
                        affected_parity_qubits.insert(q);
                    }
                }
            }
            std::cout << "\n";
        }

        const uint max_level = curr_max_detector / detectors_per_round;
        for (int32_t q : affected_parity_qubits) {
            clean_parity_qubit(q, 0, max_level, s);
        }
    }
}

void
Fleece::clean_parity_qubit(uint32_t q, uint lower_level, uint upper_level, uint64_t shot) {
    if (fake_run) {
        return;
    }

    auto v = lattice_graph.get_vertex_by_qubit(q);
    const uint8_t offset = v->base_detector < (detectors_per_round>>1) ? 0 : 1;
    if (offset && lower_level == 0) {
        return;
    }
    lower_level -= offset;
    upper_level -= offset;

    uint8_t prev_meas;
    if (lower_level == 0) {
        prev_meas = 0;
    } else {
        uint32_t prev_meas_t = v->measurement_times[lower_level - 1];
        prev_meas = sim.m_record.storage[prev_meas_t][shot];
    }
    for (uint i = lower_level; i < upper_level; i++) {
        uint32_t meas_t = v->measurement_times[i];
        sim.m_record.storage[meas_t][shot] = prev_meas;
        sim.leak_record.storage[meas_t][shot] = 0;
    }
    // Clear syndrome bits.
    uint currd = v->base_detector + lower_level*detectors_per_round;
    const uint maxd = v->base_detector + upper_level*detectors_per_round;
    while (currd < maxd) {
        meas_results[currd][shot] = 0;
        leak_results[currd][shot] = 0;
        currd = next_detector(currd);
    }
}

uint
Fleece::base_detector(uint detector) {
    return detector % detectors_per_round;
}

uint
Fleece::next_detector(uint detector) {
    return detector + detectors_per_round;
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
        }
    }
    for (uint i = 0; i < num_observables; i++) {
        mod[mod_ptr++] |= orig[orig_ptr++];
    }
    return mod;
}

} // qrc
