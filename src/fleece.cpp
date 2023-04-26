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
    n_lrus_used(0),
    n_leaks_removed_by_lru(0),
    n_leaks_removed_were_visible(0),
    parity_leakage_population(params.rounds, 0),
    data_leakage_population(params.rounds, 0),
    circuit_params(params),
    flags(flags),
    reset_basis(reset_basis),
    output_basis(output_basis),
    base_circuit(stim::generate_surface_code_circuit(params).circuit),
    sim(nullptr),
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
Fleece::create_syndromes(uint64_t shots, bool maintain_failure_log, bool record_in_rtanalyzer, uint periodicity) {
    const size_t num_results = base_circuit.count_detectors() + base_circuit.count_observables();
    stim::simd_bit_table syndromes(2*shots, num_results);
    syndromes.clear();

    const uint64_t rng_mod = (uint32_t) (1.0/(circuit_params.before_measure_flip_probability));

    uint64_t syndrome_table_index = 0;

    const uint dist = circuit_params.distance;
    const uint detectors_per_round = (dist*dist - 1) >> 1;

    if (periodicity < 2) {
        periodicity = dist;
    }
    while (shots) {
        leakage_enabled = true;
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
            if ((output_basis == 'X') == v->is_x_parity) {
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
        uint mem_depth = (dist - 1) >> 1;
        std::vector<std::vector<uint8_t>> prev_syndromes(mem_depth, std::vector<uint8_t>(parity_qubits.size(), 0));
        std::vector<std::vector<uint8_t>> prev_leakages(mem_depth, std::vector<uint8_t>(parity_qubits.size(), 0));

        std::set<fleece::LatticeGraph::Vertex*> already_swapped;

        std::deque<fleece::LatticeGraph::Vertex*> usage_queue;
        for (auto v : data_qubits) {
            usage_queue.push_back(v);
        }

        for (auto v : parity_qubits) {
            usage_queue.push_back(v);
        }

        std::set<fleece::LatticeGraph::Vertex*> involved_in_lrc;

        uint32_t meas_t_offset = parity_qubits.size();
        for (uint r = 1; r <= circuit_params.rounds; r++) {
            std::set<fleece::LatticeGraph::Vertex*> infected;   // Set of vertices with potential leakages.
            // Check for problems in the prior round.
            const uint n = parity_qubits.size();

            uint index = r - 1;
            for (auto v : data_qubits) {
                dlp[index] += sim->leakage_table[v->qubit][0];
            }
            for (auto v : parity_qubits) {
                plp[index] += sim->leakage_table[v->qubit][0];
            }
            
            if (r == circuit_params.rounds) break;

            if (r == circuit_params.rounds - 1 && (flags & LAST_R_NO_LEAKAGE)) {
                leakage_enabled = false;
                for (uint q = 0; q < base_circuit.count_qubits(); q++) {
                    if (sim->leakage_table[q][0]) {
                        sim->leakage_table[q][0] = 0;
                        sim->x_table[q][0] = rng() & 0x1;
                        sim->z_table[q][0] = rng() & 0x1;
                    }
                }
            }

            std::vector<uint8_t> comb_syndrome(syndrome);
            for (uint i = 0; i < prev_syndromes.size(); i++) {
                auto ps = prev_syndromes[i];
                auto pl = prev_leakages[i];
                for (uint i = 0; i < comb_syndrome.size(); i++) {
                    if (M_LRC_L_RESET_3WAY) {
                        comb_syndrome[i] |= ps[i] & ~pl[i];
                    } else {
                        comb_syndrome[i] |= ps[i];
                    }
                }
            }

            // Check for infections or leakages in the current round.
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                if (!(flags & (MARS_MITIGATION | M_MLR_W_ALAP_CORR))) {
                    continue;
                }
                if (!syndrome[i]) {
                    continue;
                }
                auto v = parity_qubits[i];
                uint k = circuit_params.distance >> 1;
                for (auto w : lattice_graph.adjacency_list(v)) {
                    if (w == acting_parity_qubits[stab_meas_time[v]]) {
                        uint adj_flips = 0;
                        for (auto u : lattice_graph.get_orderk_neighbors(v, 2)) {
                            if (!u->is_data && comb_syndrome[stab_meas_time[u]]) {
                                adj_flips++;
                            }
                        }
                        if (adj_flips > 2 || (adj_flips & 0x1)) {
                            infected.insert(v);
                        }
                    } else {
                        uint adj_flips = 0;
                        for (auto u : lattice_graph.adjacency_list(w)) {
                            if (u != v && comb_syndrome[stab_meas_time[u]]) {
                                adj_flips++;
                            }
                        }
                        if (adj_flips) {
                            infected.insert(w);
                        } else if (leakages[stab_meas_time[v]] && (flags & M_LRC_L_RESET_3WAY)) {
                            infected.insert(w);
                        }
                    }
                }
            }

            // Perform mitigative action (SWAP LRCs, clearing infected qubits, etc.)
            std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*> swap_targets;
            std::set<fleece::LatticeGraph::Vertex*> used;
            if ((flags & (MARS_MITIGATION | M_MLR_W_ALAP_CORR)) && !(flags & M_PERIODICITY)) {
                for (auto v : infected) {
                    if (v->is_data && swap_set.count(v)) {
                        auto w = swap_set[v];
                        if (!infected.count(w)) {
                            swap_targets[w] = v;
                        }
                    }
                }

                if (infected.count(unlucky_data_qubit)) {
                    for (auto w : lattice_graph.adjacency_list(unlucky_data_qubit)) {
                        if (swap_targets.count(w) || infected.count(w)) {
                            continue;
                        }
                        swap_targets[w] = unlucky_data_qubit;
                        break;
                    }
                }
            } else if (flags & M_PERIODICITY) {
                if (r % periodicity == periodicity - 1) {
                    // Use SWAP LRCs on unswapped qubits.
                    for (auto v : data_qubits) {
                        if (!already_swapped.count(v) && swap_set.count(v)) {
                            swap_targets[swap_set[v]] = v;
                        }
                    }
                }

                for (auto v : infected) {
                    if (v->is_data && swap_set.count(v)) {
                        auto w = swap_set[v];
                        if (!swap_targets.count(w) && !infected.count(w)) {
                            swap_targets[w] = v;
                            already_swapped.insert(v);
                        }
                    }
                }

                if (infected.count(unlucky_data_qubit)) {
                    for (auto w : lattice_graph.adjacency_list(unlucky_data_qubit)) {
                        if (!swap_targets.count(w) && !infected.count(w)) {
                            swap_targets[w] = unlucky_data_qubit;
                            already_swapped.insert(w);
                            break;
                        }
                    }
                }

                if (r % periodicity == periodicity - 1) {
                    already_swapped.clear();
                }
            } else if (flags & M_IDEAL_LRU) {
                std::vector<fleece::LatticeGraph::Vertex*> leaked_data_qubits;
                std::vector<fleece::LatticeGraph::Vertex*> unleaked_parity_qubits;
                
                for (auto v : data_qubits) {
                    if (sim->leakage_table[v->qubit][0]) {
                        leaked_data_qubits.push_back(v);
                    }
                }

                for (auto v : parity_qubits) {
                    if (!sim->leakage_table[v->qubit][0]) {
                        unleaked_parity_qubits.push_back(v);
                    }
                }

                uint nd = leaked_data_qubits.size();
                uint np = unleaked_parity_qubits.size();
                uint number_of_nodes = nd > np ? 2*nd : 2*np;

                PerfectMatching pm(number_of_nodes, 4*number_of_nodes);
                pm.options.verbose = false;

                for (uint i = 0; leaked_data_qubits.size(); i++) {
                    auto v = leaked_data_qubits[i];
                    for (uint j = 0; j < unleaked_parity_qubits.size(); j++) {
                        auto w = unleaked_parity_qubits[j];
                        for (auto u : lattice_graph.adjacency_list(v)) {
                            if (u == w) {
                                pm.AddEdge(i, nd + j, 1);
                                break;
                            }
                        }
                    }
                }

                for (uint k = nd+np; k < number_of_nodes; k++) {
                    if (nd > np) {
                        for (uint i = 0; i < leaked_data_qubits.size(); i++) {
                            pm.AddEdge(nd+i, k, 1);
                        }
                    } else {
                        for (uint i = 0; i < unleaked_parity_qubits.size(); i++) {
                            pm.AddEdge(i, k, 1);
                        }
                    }
                }

                pm.Solve();
                for (uint i = 0; i < leaked_data_qubits.size(); i++) {
                    uint jj = pm.GetMatch(i);
                    if (jj >= nd+np) {
                        // This is a dummy node.
                        i++;
                        continue;
                    }
                    auto v = leaked_data_qubits[i];
                    uint j = jj - nd;
                    auto w = unleaked_parity_qubits[j];
                    
                    swap_targets[w] = v;
                }
            } else if (flags & EN_SWAP_LRU) {
                for (auto v : usage_queue) {
                    if (used.count(v)) {
                        continue;
                    }

                    if (v->is_data) {
                        fleece::LatticeGraph::Vertex * w;
                        if (v == unlucky_data_qubit) {
                            w = lattice_graph.adjacency_list(v)[0];
                        } else {
                            w = swap_set[v];
                        }
                        if (swap_targets.count(w)) {
                            continue;
                        }
                        swap_targets[w] = v;
                    } else {
                        if (swap_targets.count(v)) {
                            continue;
                        }
                    }
                    used.insert(v);
                    if (used.size() == parity_qubits.size()) {
                        break;
                    }
                }
                std::deque<fleece::LatticeGraph::Vertex*> new_entries;
                for (auto it = usage_queue.begin(); it != usage_queue.end(); ) {
                    if (used.count(*it)) {
                        new_entries.push_back(*it);
                        used.erase(*it);
                        it = usage_queue.erase(it);
                    } else {
                        it++;
                    }
                }

                for (auto v : new_entries) {
                    usage_queue.push_back(v);
                }
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
                if (involved_in_lrc.count(v)) {
                    apply_round_start_error(v->qubit, 2);
                } else {
                    apply_round_start_error(v->qubit);
                }
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

            involved_in_lrc.clear();
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                auto v = parity_qubits[i];
                if (v->is_x_parity) {
                    apply_H(v->qubit);
                }
                
                if (swap_targets.count(v) && leakage_enabled) {
                    auto w = swap_targets[v];
                    apply_SWAP(v->qubit, w->qubit);
                    apply_measure(w->qubit);
                    apply_reset(w->qubit);

                    n_lrus_used++;
                    acting_parity_qubits[i] = w;
                } else {
                    apply_measure(v->qubit);
                    apply_reset(v->qubit);

                    acting_parity_qubits[i] = v;
                }
                uint32_t pmt = measurement_time - parity_qubits.size();
                // Model measurement error rate on leakage.
                leakages[i] = sim->leak_record.storage[measurement_time][0]; //^ ((rng() % rng_mod) == 0);
                syndrome[i] = sim->m_record.storage[measurement_time][0] ^ sim->m_record.storage[pmt][0];

                if (swap_targets.count(v) && leakages[i]) {
                    n_leaks_removed_by_lru++;
                    auto w = swap_targets[v];
                    for (auto u : lattice_graph.adjacency_list(w)) {
                        if (prev_syndromes[0][stab_meas_time[u]]) {
                            n_leaks_removed_were_visible++;
                            break;
                        }
                    }
                }

                // Apply error to 3 way measurement
                if (flags & M_LRC_L_RESET_3WAY) {
                    leakages[i] ^= ((rng() % rng_mod) == 0);
                }

                if (swap_targets.count(v)) {
                    auto w = swap_targets[v];
                    if (leakages[i] && (flags & M_LRC_L_RESET_3WAY)) {
                        // Adaptively perform a second reset when detecting a leakage error.
                        apply_reset(v->qubit);
                    } else {
                        // No need to apply error on "swap back" -- we can modify the next syndrome extraction
                        // round, but this is easier to model.
                        apply_SWAP(v->qubit, w->qubit, 2);
                        involved_in_lrc.insert(v);
                        involved_in_lrc.insert(w);
                    }
                }
                measurement_time++;
            }

            /*
            if ((flags & MARS_MITIGATION) && !(flags & M_LRC_L_RESET_3WAY)) {
                for (auto v : parity_qubits) {
                    if (!swap_targets.count(v)) {
                        leakages[stab_meas_time[v]] = 0;
                        continue;
                    }
                    auto w = swap_targets[v];
                    auto neighborhood = lattice_graph.get_orderk_neighbors(v, 2);
                    uint parity_flips = 0;
                    uint neighborhood_size = 0;
                    for (auto u : neighborhood) {
                        if (u->is_data || v == u) {
                            continue;
                        }
                        neighborhood_size++;
                        parity_flips += syndrome[stab_meas_time[u]];
                    }

                    auto adj_list = lattice_graph.adjacency_list(w);
                    uint data_flips = 0;
                    for (auto u : adj_list) {
                        if (u == v) {
                            continue;
                        }
                        data_flips += syndrome[stab_meas_time[u]] + prev_syndromes[0][stab_meas_time[u]];
                    }
                    if (parity_flips > 4 && data_flips > adj_list.size()-1) {
                        leakages[stab_meas_time[v]] = 1;
                        apply_reset(v->qubit);
                    } else {
                        leakages[stab_meas_time[v]] = 0;
                        apply_SWAP(v->qubit, w->qubit, 2);
                        involved_in_lrc.insert(v);
                        involved_in_lrc.insert(w);
                    }
                }
            }
            */

            if (flags & M_LRC_L_RESET_3WAY)  {
                for (auto v : parity_qubits) {
                    uint mt = stab_meas_time[v] + meas_t_offset;
                    uint pmt = mt - parity_qubits.size();
                    if (!leakages[stab_meas_time[v]]) {
                        continue;
                    }

                    auto neighborhood = lattice_graph.get_orderk_neighbors(v, 2);
                    uint flips = 0;
                    for (auto u : neighborhood) {
                        if (u->is_data) {
                            continue;
                        }
                        flips += (syndrome[stab_meas_time[u]] && ~leakages[stab_meas_time[u]]);
                    }

                    if (flips) {
                        syndrome[stab_meas_time[v]] = 1;
                        sim->m_record.storage[mt][0] = ~sim->m_record.storage[pmt][0];
                    } else {
                        syndrome[stab_meas_time[v]] = 0;
                        sim->m_record.storage[mt][0] = sim->m_record.storage[pmt][0];
                    }
                }
            }

            // Determine if a parity qubit has leaked during an LRC.
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
        stim::read_from_sim(*sim, det_obs, false, true, true, measure_results, leakage_results);
        measure_results = measure_results.transposed();
        leakage_results = leakage_results.transposed();
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
            swap_set[data_qubits[i]] = parity_qubits[j - data_qubits.size()];
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
    if (leakage_enabled) {
        if (roundleak_table.count(qubit)) {
            p = roundleak_table[qubit];
        } else {
            p = circuit_params.get_before_round_leakage_probability();
            roundleak_table[qubit] = p;
        }

        stim::OperationData lerror_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->LEAKAGE_ERROR(lerror_data);
    }
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

    if (leakage_enabled) {
        if (cliffordleak_table.count(q1_q2)) {
            p = cliffordleak_table[q1_q2];
        } else {
            p = circuit_params.get_after_clifford_leakage_probability();
            cliffordleak_table[q1_q2] = p;
        }
        stim::OperationData leak_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
        sim->LEAKAGE_ERROR(leak_data);
    }

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
Fleece::apply_SWAP(uint32_t qubit1, uint32_t qubit2, uint cnots, bool add_error) {
    for (uint i = 0; i < cnots; i++) {
        fp_t p = 0.0;
        stim::GateTarget q1{qubit1};
        stim::GateTarget q2{qubit2};
        std::vector<stim::GateTarget> gate_targets;
        if (i & 0x1) {
            gate_targets = std::vector<stim::GateTarget>{q2, q1};
        } else {
            gate_targets = std::vector<stim::GateTarget>{q1, q2};
        }
        stim::OperationData cnot_data{PTR<double>(&p), PTR<stim::GateTarget>(gate_targets)};
        sim->ZCX(cnot_data);

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

            if (leakage_enabled) {
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
