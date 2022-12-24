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
        std::cout << "Round = " << round << "\n";
        std::cout << "(min, max, ndet) = (" 
            << curr_min_detector << ", " << curr_max_detector 
            << ", " << circuit.count_detectors() << ")\n";
        std::cout << "leakage enabled: " << sim.leakage_enabled << "\n";
#endif
        stim::read_from_sim(sim, det_obs, false, false, true, 
                meas_results, leak_results, curr_max_detector);
        // Correct data qubit "tower-like" errors.
        if (!fake_run) {
            tower_correct(shots);
        }

        // Update curr min and curr max detectors.
        curr_min_detector += detectors_per_round; 
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
#define TOLERANCE 2

void
Fleece::tower_correct(uint64_t shots) {
    typedef std::array<uint8_t, TOWER_SIZE> tower_t;
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
        std::map<int32_t, tower_t> tower_table;
        for (uint i = 0; i < curr_max_detector; i++) {
            if (!meas_results_T[s][i]) {
                continue;
            }
            uint level = get_level(i);
            auto pqv = lattice_graph.get_vertex_by_detector(base_detector(i));
            for (auto w : lattice_graph.adjacency_list(pqv)) {
                if (!tower_table.count(w->qubit)) {
                    tower_table[w->qubit] = tower_t();
                    tower_table[w->qubit].fill(0);
                }
                tower_table[w->qubit][level]++;
            }
        }

        const uint32_t max_level = curr_max_detector / detectors_per_round;
        for (auto pair : tower_table) {
            int32_t q = pair.first;
            tower_t det_tower = pair.second;

            auto dqv = lattice_graph.get_vertex_by_qubit(q);
            auto adjlist = lattice_graph.adjacency_list(dqv);
            const bool is_next_to_boundary = adjlist.size() < 4;
            
            uint tower_weight = 0;
            uint8_t zero_levels = 0;
            int32_t lower_level = max_level;
            while (lower_level >= 0) {
                // If the qubit is next to the boundary, then we expect to only
                // flip one parity qubit on a leakage 50% of the time and
                // both parity qubits 25% of the time.
                //
                // Otherwise, we expect to flip 2 parity qubits 50% of the time
                // and 4 parity qubits 25% of the time.
                if (det_tower[lower_level] >= (1+!is_next_to_boundary)) {
                    tower_weight++;
                    zero_levels = 0;
                } else {
                    zero_levels++;
                }
                if (zero_levels > TOLERANCE) {
                    break;
                }
                lower_level--;
            }

            if (lower_level < 0) {
                lower_level = 0;
            }
            
            if (tower_weight >= tower_cutoff) {
#ifdef FLEECE_DEBUG
                std::cout << "\t[ shot " << s << " ]"
                            << " Corrected detected leakage on "
                            << "qubit " << q << "\n";
                std::cout << "\t\tis_leaked = " << sim.leakage_table[q][s] << "\n";
                std::cout << "\t\ttower_cutoff = " <<  tower_cutoff 
                            << ", tower weight = " << tower_weight << "\n";
#endif
                // Reset the data qubit.
                sim.x_table[q][s] = 0;
                if (sim.guarantee_anticommutation_via_frame_randomization) {
                    sim.z_table[q][s] = LOCAL_RNG() & 0x1;
                }
                sim.leakage_table[q][s] = 0;
                for (auto w : adjlist) {
                    if (w->qubit != LATTICE_BOUNDARY) {
                        clean_parity_qubit(w->qubit, lower_level, max_level, s);
                    }
                }
            }
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
    for (uint i = lower_level; i <= upper_level; i++) {
        uint32_t meas_t = v->measurement_times[i];
        sim.m_record.storage[meas_t][shot] = prev_meas;
        sim.leak_record.storage[meas_t][shot] = 0;
    }
    // Clear syndrome bits.
    uint currd = v->base_detector + lower_level*detectors_per_round;
    const uint maxd = v->base_detector + upper_level*detectors_per_round;
    while (currd <= maxd) {
        meas_results[currd][shot] = 0;
        leak_results[currd][shot] = 0;
        currd = next_detector(currd);
    }
}

uint
Fleece::get_level(uint detector) {
    return (detector + (detectors_per_round>>1)) / detectors_per_round;
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
