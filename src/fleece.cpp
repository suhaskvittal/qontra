/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#include "fleece.h"

namespace qrc {

template <typename T>
using PTR = stim::ConstPointerRange<T>;

Fleece::Fleece(const stim::CircuitGenParameters& params,
                uint16_t flags,
                std::mt19937_64& rng,
                char reset_basis,
                char output_basis)
    :failure_log(),
    rtanalyzer(nullptr),
    n_restarts(0),
    parity_leakage_population(params.rounds, 0),
    data_leakage_population(params.rounds, 0),
    circuit_params(params),
    flags(flags),
    reset_basis(reset_basis),
    output_basis(output_basis),
    base_circuit(stim::generate_surface_code_circuit(params).circuit),
    sim(nullptr),
    decoding_graph(),
    decoder_path_table(),
    lattice_graph(),
    data_qubits(),
    swap_set(),
    unlucky_data_qubit(),
    parity_qubits(),
    rounddp_table(),
    clifforddp1_table(),
    clifforddp2_table(),
    premeasflip_table(),
    postresflip_table(),
    roundleak_table(),
    cliffordleak_table(),
    leaktransport_table(),
    rng(rng)
{   
    decoding_graph = to_decoding_graph(base_circuit);
    decoder_path_table = compute_path_table(decoding_graph); 
    rtanalyzer = new fleece::RealTimeAnalyzer(base_circuit, rng);

    sim = new stim::FrameSimulator(base_circuit.count_qubits(), 1, SIZE_MAX, rng);
    lattice_graph = fleece::to_lattice_graph(base_circuit);
    lattice_path_table = fleece::compute_path_table(lattice_graph);
    
    std::vector<fleece::LatticeGraph::Vertex*> x_parity_qubits;
    std::vector<fleece::LatticeGraph::Vertex*> z_parity_qubits;
    for (auto v : lattice_graph.vertices()) {
        if (v->qubit < 0) {
            continue;
        }

        if (v->is_data) {
            data_qubits.push_back(v);
        } else if (v->is_x_parity) {
            x_parity_qubits.push_back(v);
        } else {
            z_parity_qubits.push_back(v);
        }
    }

    struct _vertex_cmp {
        typedef fleece::LatticeGraph::Vertex V;
        bool operator()(V * v1, V * v2) const {
            return v1->qubit < v2->qubit;
        }
    } cmp;

    std::sort(data_qubits.begin(), data_qubits.end(), cmp);
    std::sort(x_parity_qubits.begin(), x_parity_qubits.end(), cmp);
    std::sort(z_parity_qubits.begin(), z_parity_qubits.end(), cmp);

    for (auto v : x_parity_qubits) {
        parity_qubits.push_back(v);
    }

    for (auto v : z_parity_qubits) {
        parity_qubits.push_back(v);
    }

    compute_optimal_swap_set();
}

Fleece::~Fleece() {
    delete rtanalyzer;
    delete sim;
}

stim::simd_bit_table
Fleece::create_syndromes(uint64_t shots, bool maintain_failure_log, bool record_in_rtanalyzer) {
    const size_t num_results = base_circuit.count_detectors() + base_circuit.count_observables();
    stim::simd_bit_table syndromes(2*shots, num_results);
    syndromes.clear();

    const uint64_t rng_mod = (uint32_t) (1.0/circuit_params.before_measure_flip_probability);

    uint64_t syndrome_table_index = 0;
    bool restart_shot = false;
    while (shots) {
        std::vector<uint64_t> dlp(data_leakage_population.size(), 0);
        std::vector<uint64_t> plp(parity_leakage_population.size(), 0);

        std::string shot_log;
        sim->reset_all();
        // Do first round -- special.
        //
        // Reset all qubits.
        uint32_t measurement_time = 0;
        std::vector<uint8_t> syndrome(parity_qubits.size(), 0);
        std::vector<uint8_t> leakages(parity_qubits.size(), 0);
        std::vector<fleece::LatticeGraph::Vertex*> acting_parity_qubits(parity_qubits.size());

        std::map<fleece::LatticeGraph::Vertex*, uint32_t> stab_meas_time;
            
        for (auto v : data_qubits) {
            apply_reset(v->qubit);
            if (output_basis == 'X') {
                apply_H(v->qubit, false);
            }
            apply_round_start_error(v->qubit);
        }
        for (auto v : parity_qubits) {
            apply_reset(v->qubit);
            if (v->is_x_parity) {
                apply_H(v->qubit);
            }
        }
        // Perform stabilizer measurements.
        for (uint j = 0; j < 4; j++) {
            for (auto v : parity_qubits) {
                auto w = lattice_graph.get_cx_mate(v, j);
                if (w == nullptr) {
                    continue;
                }

                if (v->is_x_parity) {
                    apply_CX(v->qubit, w->qubit);
                } else {
                    apply_CX(w->qubit, v->qubit);
                }
            }
        }
        
        for (auto v : parity_qubits) {
            if (v->is_x_parity) {
                apply_H(v->qubit);
            }
            apply_measure(v->qubit);
            apply_reset(v->qubit);

            leakages[measurement_time] = sim->leak_record.storage[measurement_time][0] ^ ((rng() % rng_mod) == 0);
            if (v->is_x_parity) {
                syndrome[measurement_time] = 0;
            } else {
                // Only record as error if parity qubit isn't leaked.
                syndrome[measurement_time] = sim->m_record.storage[measurement_time][0];
            }
            acting_parity_qubits[measurement_time] = v;
            stab_meas_time[v] = measurement_time;
            measurement_time++;
        }
        // Write syndromes to failure log.
        if (maintain_failure_log) {
            shot_log += "R0------------------------\n";
            write_leakage_condition_to_log(shot_log);
            shot_log += "\taliases:";
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                shot_log += "\tP" + std::to_string(parity_qubits[i]->qubit);
            }
            shot_log += "\n";
            write_syndrome_to_log(shot_log, syndrome, leakages);
        }
        
        if (record_in_rtanalyzer) {
            for (auto v : lattice_graph.vertices()) {
                if (v->qubit < 0) {
                    continue;
                }

                if (sim->leakage_table[v->qubit][0]) {
                    rtanalyzer->increment(v, syndrome_table_index>>1);
                } else {
                    rtanalyzer->read_and_flush(v, syndrome_table_index>>1);
                }
            }
        }
        // Now, we need to simulate the circuit operation by operation, round by round.
        std::array<std::vector<uint8_t>, 2> prev_syndromes;
        std::array<std::vector<uint8_t>, 2> prev_leakages;

        std::vector<uint> matched_detectors;  // Already matched and fixed -- update syndrome at the end.
        std::vector<uint> result_obs(base_circuit.count_observables(), 0);  // Observable resulting from matchings.

        prev_syndromes.fill(std::vector<uint8_t>(parity_qubits.size(), 0));
        prev_leakages.fill(std::vector<uint8_t>(parity_qubits.size(), 0));

        std::set<fleece::LatticeGraph::Vertex*> already_swapped;

        std::deque<fleece::LatticeGraph::Vertex*> usage_queue;
        for (auto v : data_qubits) {
            usage_queue.push_back(v);
        }
        for (auto v : parity_qubits) {
            usage_queue.push_back(v);
        }

        std::deque<fleece::LatticeGraph::Vertex*> decoding_queue;

        uint32_t meas_t_offset = 0;
        restart_shot = false;
        for (uint r = 1; r <= circuit_params.rounds; r++) {
            std::set<fleece::LatticeGraph::Vertex*> infected;   // Set of vertices with potential leakages.
            // Check for problems in the prior round.
            std::set<fleece::LatticeGraph::Vertex*> faulty_parity_checks;
            const uint n = parity_qubits.size();

            uint index = r - 1;
            for (auto v : data_qubits) {
                dlp[index] += sim->leakage_table[v->qubit][0];
            }
            for (auto v : parity_qubits) {
                plp[index] += sim->leakage_table[v->qubit][0];
            }

            // Dump decoding queue.
            while (!decoding_queue.empty()) {
                auto v = decoding_queue.front();
                decoding_queue.pop_front();

                std::vector<fleece::LatticeGraph::Vertex*> major_faults;
                if (v->is_data) {
                    major_faults.push_back(v);
                } else {
                    major_faults = lattice_graph.adjacency_list(v);
                }

                const uint dist = circuit_params.distance;
                const uint detectors_per_round = (dist*dist - 1) >> 1;

                const uint syndrome_size = base_circuit.count_detectors() + base_circuit.count_observables();
                std::vector<uint8_t> micro_syndrome(syndrome_size, 0);

                uint min_r = r > prev_syndromes.size() ? r - prev_syndromes.size() : 0;
                uint k = min_r * detectors_per_round;
    
                for (uint d : matched_detectors) {
                    if (d >= k) {
                        micro_syndrome[d] ^= 1;
                    }
                }

                for (uint syn_index = 0; syn_index < r - min_r; syn_index++) { 
                    auto syn = prev_syndromes[prev_syndromes.size() - syn_index - 1];
                    auto lk = prev_leakages[prev_leakages.size() - syn_index - 1];
                    for (uint i = 0; i < n; i++) {
                        if ((parity_qubits[i]->is_x_parity && output_basis == 'Z')
                            || (!parity_qubits[i]->is_x_parity && output_basis == 'X')) 
                        {
                            continue;
                        }
                        micro_syndrome[k++] ^= syn[i] & ~lk[i];
                    }
                }

                for (uint i = 0; i < n; i++) {
                    if ((parity_qubits[i]->is_x_parity && output_basis == 'Z')
                        || (!parity_qubits[i]->is_x_parity && output_basis == 'X')) 
                    {
                        continue;
                    }
                    micro_syndrome[k++] ^= syndrome[i] & ~leakages[i];
                }

                auto mdres = decode_error(micro_syndrome, major_faults);
                for (auto pair : mdres.matching) {
                    if (pair.first > pair.second) {
                        continue;
                    }
                    matched_detectors.push_back(pair.first);
                    if (pair.second != BOUNDARY_INDEX) {
                        matched_detectors.push_back(pair.second);
                    }
                }

                for (uint i = 0; i < base_circuit.count_observables(); i++) {
                    result_obs[i] ^= mdres.correction[i];
                }
            }

            if (r == circuit_params.rounds) break;

            // Check for infections or leakages in the current round.
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                if (flags & NO_MITIGATION) continue;
                if (leakages[i]) {
                    auto v = acting_parity_qubits[i];
                    decoding_queue.push_back(v);
                } else if (syndrome[i]) {
                    auto v = parity_qubits[i];
                    for (uint32_t j = i+1; j < parity_qubits.size(); j++) {
                        if (syndrome[j]) {
                            auto w = parity_qubits[j];
                            auto path = lattice_path_table[std::make_pair(v, w)].path;
                            uint chain_length = path.size() >> 1;
                            if (chain_length < ((circuit_params.distance-1) >> 1)) {
                                for (auto u : path) infected.insert(u);
                            }
                        }
                    }
                }
            }
            // Using infected set, identify SWAP LRCs.
            std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*> swap_targets;
            if ((r % circuit_params.distance == circuit_params.distance - 1) && (flags & EN_STATE_DRAIN)) {
                // Perform last minute swap LRUs to kill any remaining leakage errors.
                for (auto pair : swap_set) {
                    swap_targets[pair.first] = pair.second;
                }
                // We will simulate "Temporary Stabilizer Extension" on the unlucky data qubit.
            } else if (flags & EN_SWAP_LRU) {
                if (r % 3 == 1) {
                    already_swapped.clear();
                }
                
                for (auto w : parity_qubits) {
                    if (swap_targets.count(w)) {
                        continue;
                    }
                    for (auto v : lattice_graph.adjacency_list(w)) {
                        if (already_swapped.count(v)) {
                            continue;
                        }
                        swap_targets[w] = v;
                        already_swapped.insert(v);
                    }
                }
            } else if (!(flags & NO_MITIGATION)) {
                std::map<fleece::LatticeGraph::Vertex*, uint> priority_table;
                for (uint i = 0; i < usage_queue.size(); i++) {
                    priority_table[usage_queue[i]] = i;
                }
                std::set<fleece::LatticeGraph::Vertex*> swapped;
                for (auto v : infected) {
                    if (already_swapped.count(v)) {
                        continue;
                    }
                    int max_priority = -1;
                    fleece::LatticeGraph::Vertex * victim;
                    for (auto w : lattice_graph.adjacency_list(v)) {
                        if (swap_targets.count(w)) {
                            continue;
                        }
                        uint prio = priority_table[w];
                        if (prio > max_priority) {
                            max_priority = prio;
                            victim = w;
                        }
                    }
                    swap_targets[victim] = v;
                    swapped.insert(v);
                }

                for (auto v : usage_queue) {
                    if (swap_targets.count(v) || swapped.count(v)) {
                        continue;
                    }
                    swapped.insert(v);
                    if (swapped.size() == parity_qubits.size()) {
                        break;
                    }
                }

                std::deque<fleece::LatticeGraph::Vertex*> new_queue_entries;
                for (auto it = usage_queue.begin(); it != usage_queue.end(); ) {
                    if (swapped.count(*it)) {
                        it = usage_queue.erase(it);
                        new_queue_entries.push_back(*it);
                    } else {
                        it++;
                    }
                }

                for (auto v : new_queue_entries) {
                    usage_queue.push_back(v);
                }
                already_swapped = swapped;
            }

            // Update previous syndromes and leakages
            for (uint j = prev_syndromes.size() - 1; j >= 1; j--) {
                prev_syndromes[j] = prev_syndromes[j-1];
                prev_leakages[j] = prev_leakages[j-1];
            }
            prev_syndromes[0] = syndrome;
            prev_leakages[0] = leakages;

            // Perform operations.
            for (auto v : data_qubits) {
                apply_round_start_error(v->qubit);
            }

            for (auto v : parity_qubits) {
                if (v->is_x_parity) {
                    apply_H(v->qubit);
                }
            }
            // Perform stabilizer measurements.
            for (uint j = 0; j < 4; j++) {
                for (auto v : parity_qubits) {
                    auto w = lattice_graph.get_cx_mate(v, j);
                    if (w == nullptr) {
                        continue;
                    }

                    if (v->is_x_parity) {
                        apply_CX(v->qubit, w->qubit);
                    } else {
                        apply_CX(w->qubit, v->qubit);
                    }
                }
            }

            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                auto v = parity_qubits[i];
                if (v->is_x_parity) {
                    apply_H(v->qubit);
                }
                
                if (swap_targets.count(v)) {
                    auto _v = swap_targets[v];
                    apply_SWAP(v->qubit, _v->qubit);
                    apply_measure(_v->qubit);
                    apply_reset(_v->qubit);
                    // Don't apply a depolarizing error if this is the last syndrome extraction round.
                    // Technically, we could just measure the parity qubits to measure the data qubits,
                    // but this is a bit harder.
                    apply_SWAP(v->qubit, _v->qubit, (r != circuit_params.rounds - 1));

                    acting_parity_qubits[i] = _v;
                } else {
                    apply_measure(v->qubit);
                    apply_reset(v->qubit);

                    acting_parity_qubits[i] = v;
                }
                uint32_t pmt = measurement_time - parity_qubits.size();
                // Update syndromes.
                // Model measurement error rate on leakage.
                leakages[i] = sim->leak_record.storage[measurement_time][0] ^ ((rng() % rng_mod) == 0);
                if (leakages[i] && swap_targets.count(v)) {
                    if (r == circuit_params.rounds-1) {
                        restart_shot = true;
                    }
                }
                syndrome[i] = sim->m_record.storage[measurement_time][0] ^ sim->m_record.storage[pmt][0];

                measurement_time++;
            }

            if ((flags & EN_TMP_STAB_EXT)) {
                // Simulate "Temporary Stabilizer Extension"
                // We apply the same errors as we would at the beginning of a round (depolarizing + leakage)
                apply_round_start_error(unlucky_data_qubit->qubit);
                if (sim->leakage_table[unlucky_data_qubit->qubit][0] ^ ((rng() % rng_mod) == 0)) {
                    restart_shot = (r == circuit_params.rounds - 1);
                }
                if (sim->leakage_table[unlucky_data_qubit->qubit][0]) {
                    // Only situation where the SWAP fails.
                    apply_reset(unlucky_data_qubit->qubit);
                }
            }

            if (maintain_failure_log) {
                shot_log += "R" + std::to_string(r) + "------------------------\n";
                write_leakage_condition_to_log(shot_log);
                write_aliases_to_log(shot_log, swap_targets);
                write_syndrome_to_log(shot_log, syndrome, leakages);
            }

            if (record_in_rtanalyzer) {
                for (auto v : lattice_graph.vertices()) {
                    if (v->qubit < 0) {
                        continue;
                    }

                    if (sim->leakage_table[v->qubit][0]) {
                        rtanalyzer->increment(v, syndrome_table_index>>1);
                    } else {
                        rtanalyzer->read_and_flush(v, syndrome_table_index>>1);
                    }
                }
            }

            meas_t_offset += parity_qubits.size();
        }
        // Perform tail measurements.
        if (restart_shot && !(flags & (NO_MITIGATION | NO_RESTART))) {
            n_restarts++;
            continue;
        }

        for (uint i = 0; i < data_leakage_population.size(); i++) {
            data_leakage_population[i] += dlp[i];
            parity_leakage_population[i] += plp[i];
        }

        for (auto v : data_qubits) {
            if (output_basis == 'X') {
                apply_H(v->qubit, false);
            }
            apply_measure(v->qubit);
        }
        // Everything is done, extract the syndrome.
        stim::simd_bit_table measure_results(num_results, 1);
        stim::simd_bit_table leakage_results(num_results, 1);
        auto det_obs = stim::DetectorsAndObservables(base_circuit);
        stim::read_from_sim(*sim, det_obs, false, true, true,
                                measure_results, leakage_results);
        measure_results = measure_results.transposed();
        leakage_results = leakage_results.transposed();
        for (uint d : matched_detectors) {
            measure_results[0][d] ^= 1;
        }
        for (uint i = 0; i < result_obs.size(); i++) {
            measure_results[0][i+base_circuit.count_detectors()] ^= result_obs[i];
        }
        syndromes[syndrome_table_index++] |= measure_results[0];
        syndromes[syndrome_table_index++] |= leakage_results[0];
        shots--;

        if (leakage_results[0][base_circuit.count_detectors()] && maintain_failure_log) {
            failure_log += "===============================\n" + shot_log;
        }
    }

    if (record_in_rtanalyzer) {
        rtanalyzer->flush_table();
    }
    return syndromes;
}

DecoderShotResult
Fleece::decode_error(const std::vector<uint8_t>& syndrome, const std::vector<fleece::LatticeGraph::Vertex*>& faults) {
    std::vector<uint> detector_list;
    for (uint d = 0; d < base_circuit.count_detectors(); d++) {
        if (syndrome[d]) {
            detector_list.push_back(d);
        }
    }

    if (detector_list.size() & 0x1) {
        detector_list.push_back(BOUNDARY_INDEX);
    }

    const uint number_of_nodes = detector_list.size();
    const uint number_of_edges = (number_of_nodes*(number_of_nodes-1)) >> 1;
        
    PerfectMatching pm(number_of_nodes, number_of_edges);
    pm.options.verbose = false;
    std::set<std::pair<uint, uint>> faulting_edges;
    for (uint i = 0; i < detector_list.size(); i++) {
        uint di = detector_list[i];
        auto vi = decoding_graph.get_vertex(di);
        for (uint j = i+1; j < detector_list.size(); j++) {
            uint dj = detector_list[j];
            auto vj = decoding_graph.get_vertex(dj);

            auto w = decoder_path_table[std::make_pair(vi, vj)].distance;
            if (path_contains_faults(di, dj, faults)) {
                w -= log10(circuit_params.before_round_data_depolarization) / log10(0.5);
                if (w < 1) { 
                    w = 1.0;
                }
                faulting_edges.insert(std::make_pair(i, j));
                faulting_edges.insert(std::make_pair(j, i));
            }
            uint32_t m_weight = (uint32_t) (w * MWPM_INTEGER_SCALE);
            pm.AddEdge(i, j, m_weight);
        }
    }

    pm.Solve();

    std::map<uint, uint> matching;
    std::vector<uint8_t> correction(base_circuit.count_observables());
    for (uint i = 0; i < number_of_nodes; i++) {
        uint j = pm.GetMatch(i);
        auto i_j = std::make_pair(i, j);
        if (faulting_edges.count(i_j)) {
            matching[detector_list[i]] = detector_list[j];

            if (i < j) {
                auto vi = decoding_graph.get_vertex(detector_list[i]);
                auto vj = decoding_graph.get_vertex(detector_list[j]);

                auto path = decoder_path_table[std::make_pair(vi, vj)].path;
                for (uint i = 1; i < path.size(); i++) {
                    auto wi = path[i-1];
                    auto wj = path[i];
                    auto edge = decoding_graph.get_edge(wi, wj);
                    for (auto obs : edge->frames) {
                        if (obs >= 0) {
                            correction[obs] ^= 1;
                        }
                    }
                }
            }
        }
    }
    DecoderShotResult res = {0, 0, 0, correction, matching};
    return res;
}

bool
Fleece::path_contains_faults(uint d1, uint d2, const std::vector<fleece::LatticeGraph::Vertex*>& faults) {
    auto vd1 = decoding_graph.get_vertex(d1);
    auto vd2 = decoding_graph.get_vertex(d2);
    auto path = decoder_path_table[std::make_pair(vd1, vd2)].path;

    const uint dist = circuit_params.distance;
    const uint detectors_per_round = (dist*dist-1) >> 1;

    std::set<fleece::LatticeGraph::Vertex*> data_qubits_in_path;

    for (uint i = 1; i < path.size() - 1; i++) {
        auto wdi = path[i-1];
        auto wdj = path[i];

        if (wdi->detector == BOUNDARY_INDEX || wdj->detector == BOUNDARY_INDEX) {
            continue;
        }

        auto wli = lattice_graph.get_vertex_by_detector(wdi->detector % detectors_per_round);
        auto wlj = lattice_graph.get_vertex_by_detector(wdj->detector % detectors_per_round);
        auto common = lattice_graph.get_common_neighbors(wli, wlj);
        for (auto u : common) {
            data_qubits_in_path.insert(u);
        }
    } 

    bool has_fault = true;
    for (auto v : faults) {
        has_fault &= data_qubits_in_path.count(v);
    }
    return has_fault;
}

void
Fleece::compute_optimal_swap_set() {
    // Just get a perfect matching between data qubits and {parity qubits U dummy}
    const uint number_of_nodes = data_qubits.size() + parity_qubits.size() + 1;
    const uint dummy_number = data_qubits.size() + parity_qubits.size();
    PerfectMatching pm(number_of_nodes, number_of_nodes*4);
    pm.options.verbose = false;

    for (uint i = 0; i < data_qubits.size(); i++) {
        auto adj = lattice_graph.adjacency_list(data_qubits[i]);
        for (uint j = 0; j < parity_qubits.size(); j++) {
            // There's really no better way to check if parity_qubits[j] is in adj.
            for (auto v : adj) {
                if (parity_qubits[j] == v) {
                    pm.AddEdge(i, data_qubits.size() + j, 1);
                    break;
                }
            }
        }
        if (adj.size() == 2) {
            pm.AddEdge(i, dummy_number, 1);
        }
    }

    pm.Solve();

    for (uint i = 0; i < data_qubits.size(); i++) {
        uint j = pm.GetMatch(i);
        if (j == dummy_number) {
            unlucky_data_qubit = data_qubits[i];
        } else {
            swap_set[parity_qubits[j - data_qubits.size()]] = data_qubits[i];
        }
    }
}

void
Fleece::write_leakage_condition_to_log(std::string& shot_log) {
    shot_log += "\tLeaked qubits:";
    for (auto v : data_qubits) {
        if (sim->leakage_table[v->qubit][0]) {
            shot_log += "\tD" + std::to_string(v->qubit);
        }
    }

    for (auto v : parity_qubits) {
        if (sim->leakage_table[v->qubit][0]) {
            shot_log += "\tP" + std::to_string(v->qubit);
        }
    }
    shot_log += "\n";
}

void
Fleece::write_aliases_to_log(std::string& shot_log, 
        const std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*>& swap_targets)
{
    shot_log += "\taliases:";
    for (auto v : parity_qubits) {
        if (!swap_targets.count(v)) {
            shot_log += "\tP" + std::to_string(v->qubit);
        } else {
            const int32_t q = swap_targets.at(v)->qubit;
            shot_log += "\tD" + std::to_string(q);
        }
    }
    shot_log += "\n";
}

void
Fleece::write_syndrome_to_log(std::string& shot_log, const std::vector<uint8_t>& syndrome, 
        const std::vector<uint8_t>& leakages)
{
    shot_log += "\tsyndrome:";
    for (uint32_t i = 0; i < parity_qubits.size(); i++) {
        if (leakages[i]) {
            shot_log += "\tL";
        } else if (syndrome[i]) {
            shot_log += "\tX";
        } else {
            shot_log += "\t_";
        }
    }
    shot_log += "\n";
}

void
Fleece::apply_reset(uint32_t qubit, bool add_error) {
    fp_t p = 0.0;
    stim::GateTarget q{qubit};
    stim::OperationData reset_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->reset_z(reset_data);

    if (add_error) {
        // Apply reset error.
        if (postresflip_table.count(qubit)) {
            p = postresflip_table[qubit];
        } else {
            p = circuit_params.get_after_reset_flip_probability();
            postresflip_table[qubit] = p;
        }

        stim::OperationData xerror_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->X_ERROR(xerror_data);
    }
}

void
Fleece::apply_round_start_error(uint32_t qubit, fp_t dp_error_mult) {
    stim::GateTarget q{qubit};

    fp_t p;
    if (rounddp_table.count(qubit)) {
        p = rounddp_table[qubit];
    } else {
        p = circuit_params.get_before_round_data_depolarization();
        rounddp_table[qubit] = p;
    }
    p *= dp_error_mult;
    stim::OperationData dp_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->DEPOLARIZE1(dp_data);
    if (roundleak_table.count(qubit)) {
        p = roundleak_table[qubit];
    } else {
        p = circuit_params.get_before_round_leakage_probability();
        roundleak_table[qubit] = p;
    }

    stim::OperationData lerror_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->LEAKAGE_ERROR(lerror_data);
}

void
Fleece::apply_H(uint32_t qubit, bool add_error) {
    fp_t p = 0;
    stim::GateTarget q{qubit};
    stim::OperationData h_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->H_XZ(h_data);

    if (add_error) {
        if (clifforddp1_table.count(qubit)) {
            p = clifforddp1_table[qubit];
        } else {
            p = circuit_params.get_after_clifford_depolarization(true);
            clifforddp1_table[qubit] = p;
        }
        stim::OperationData dp_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->DEPOLARIZE1(dp_data);
    }
}

void
Fleece::apply_CX(uint32_t qubit1, uint32_t qubit2) {
    fp_t p = 0.0;
    stim::GateTarget q1{qubit1};
    stim::GateTarget q2{qubit2};
    std::vector<stim::GateTarget> gate_targets{q1, q2};
    stim::OperationData cx_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
    sim->ZCX(cx_data);

    std::vector<uint32_t> q1_q2{qubit1, qubit2};
    if (clifforddp2_table.count(q1_q2)) {
        p = clifforddp2_table[q1_q2];
    } else {
        p = circuit_params.get_after_clifford_depolarization();
        clifforddp2_table[q1_q2] = p;
    }
    stim::OperationData dp2_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
    sim->DEPOLARIZE2(dp2_data);

    if (cliffordleak_table.count(q1_q2)) {
        p = cliffordleak_table[q1_q2];
    } else {
        p = circuit_params.get_after_clifford_leakage_probability();
        cliffordleak_table[q1_q2] = p;
    }
    stim::OperationData leak_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
    sim->LEAKAGE_ERROR(leak_data);

    if (leaktransport_table.count(q1_q2)) {
        p = leaktransport_table[q1_q2];
    } else {
        p = circuit_params.get_after_clifford_leakage_transport();
    }
    stim::OperationData transport_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
    sim->LEAKAGE_TRANSPORT(transport_data);
}

void
Fleece::apply_measure(uint32_t qubit, bool add_error) {
    stim::GateTarget q{qubit};
    fp_t p;
    if (add_error) {
        if (premeasflip_table.count(qubit)) {
            p = premeasflip_table[qubit];
        } else {
            p = circuit_params.get_before_measure_flip_probability();
            premeasflip_table[qubit] = p;
        }
        stim::OperationData flip_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->X_ERROR(flip_data);
    }

    p = 0.0;
    stim::OperationData meas_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->measure_z(meas_data);
}

void
Fleece::apply_SWAP(uint32_t qubit1, uint32_t qubit2, bool add_error) {
    fp_t p = 0.0;
    stim::GateTarget q1{qubit1};
    stim::GateTarget q2{qubit2};
    std::vector<stim::GateTarget> gate_targets{q1, q2};
    stim::OperationData swap_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
    sim->SWAP(swap_data);

    if (add_error) {
        std::vector<uint32_t> q1_q2{qubit1, qubit2};
        if (clifforddp2_table.count(q1_q2)) {
            p = clifforddp2_table[q1_q2];
        } else {
            p = circuit_params.get_after_clifford_depolarization();
            clifforddp2_table[q1_q2] = p;
        }
        stim::OperationData dp2_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
        sim->DEPOLARIZE2(dp2_data);

        if (cliffordleak_table.count(q1_q2)) {
            p = cliffordleak_table[q1_q2];
        } else {
            p = circuit_params.get_after_clifford_leakage_probability();
        }
        stim::OperationData leak_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
        sim->LEAKAGE_ERROR(leak_data);

        if (leaktransport_table.count(q1_q2)) {
            p = leaktransport_table[q1_q2];
        } else {
            p = circuit_params.get_after_clifford_leakage_transport();
        }
        stim::OperationData transport_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
        sim->LEAKAGE_TRANSPORT(transport_data);
    }
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
