/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#include "fleece.h"

namespace qrc {

template <typename T>
using PTR = stim::ConstPointerRange<T>;

Fleece::Fleece(const stim::CircuitGenParameters& params,
                uint8_t flags,
                std::mt19937_64& rng,
                char reset_basis,
                char output_basis,
                bool perform_swaps)
    :circuit_params(params),
    flags(flags),
    reset_basis(reset_basis),
    output_basis(output_basis),
    perform_swaps(perform_swaps),
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
    postresleak_table(),
    rng(rng)
{
    sim = new stim::FrameSimulator(base_circuit.count_qubits(), 1, SIZE_MAX, rng);
    lattice_graph = fleece::to_lattice_graph(base_circuit);
    
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
    delete sim;
}

stim::simd_bit_table
Fleece::create_syndromes(uint64_t shots, uint disable_leakage_at_round, bool maintain_failure_log) {
    const size_t num_results = base_circuit.count_detectors() + base_circuit.count_observables();
    stim::simd_bit_table syndromes(2*shots, num_results);
    syndromes.clear();

    const uint32_t rng_mod = (uint32_t) (1.0/circuit_params.before_measure_flip_probability);

    uint64_t syndrome_table_index = 0;
    while (shots) {
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
                syndrome[measurement_time] = sim->m_record.storage[measurement_time][0] & ~leakages[measurement_time];
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
        // Now, we need to simulate the circuit operation by operation, round by round.
        std::array<std::vector<uint8_t>, 1> prev_syndrome;
        std::array<std::vector<uint8_t>, 1> prev_leakages;

        prev_syndrome.fill(std::vector<uint8_t>(parity_qubits.size(), 0));
        prev_leakages.fill(std::vector<uint8_t>(parity_qubits.size(), 0));

        uint32_t meas_t_offset = 0;
        for (uint r = 1; r < circuit_params.rounds; r++) {
            std::vector<uint8_t> comb_syndrome(parity_qubits.size(), 0);
            for (uint i = 0; i < comb_syndrome.size(); i++) {
                comb_syndrome[i] |= syndrome[i];
                for (uint j = 0; j < prev_syndrome.size(); j++) {
                    comb_syndrome[i] |= prev_syndrome[j][i];
                }
            }
            std::set<fleece::LatticeGraph::Vertex*> infected;   // Set of vertices with potential leakages.
            // Check for problems in the prior round.
            std::set<fleece::LatticeGraph::Vertex*> faulty_data_qubits_on_x;
            std::set<fleece::LatticeGraph::Vertex*> faulty_data_qubits_on_z;
            std::set<fleece::LatticeGraph::Vertex*> faulty_parity_checks;
            const uint n = parity_qubits.size();
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                if (leakages[i]) {
                    auto v = acting_parity_qubits[i];
                    // If v is a data qubit, then only the adjacent stabilizers were affected.
                    // So, clear the measurements.
                    //
                    // If v is a parity qubit, then we must also match flipped stabilizers to
                    // this leaked stabilizer measurement.
                    //
                    // TODO: account for measurement/CNOT errors on flipped stabilizers.
                    if (v->is_data) {
                        for (auto w : lattice_graph.adjacency_list(v)) {
                            faulty_parity_checks.insert(w);
                        }
                    } else if (!v->is_data) {
                        for (auto w : parity_qubits) {
                            uint32_t mt = stab_meas_time[w];
                            if (v == w || !comb_syndrome[mt]) {
                                continue;
                            }
                            auto common = lattice_graph.get_common_neighbors(v, w);
                            for (auto u : common) {
                                if (w->is_x_parity) {
                                    faulty_data_qubits_on_z.insert(u);
                                } else {
                                    faulty_data_qubits_on_x.insert(u);
                                }
                            }

                            if (!common.empty()) {
                                faulty_parity_checks.insert(w);
                            }
                        }
                    }

                } else if (syndrome[i]) {
                    // Add all data qubits involved in stabilizer measurement.
                    auto v = parity_qubits[i];
                    for (auto w : lattice_graph.adjacency_list(v)) {
                        infected.insert(w);
                    }
                }
            }

            for (auto w : faulty_data_qubits_on_x) {
                sim->x_table[w->qubit][0] ^= 1;
            }
            
            for (auto w : faulty_data_qubits_on_z) {
                sim->z_table[w->qubit][0] ^= 1;
            }

            for (auto w : faulty_parity_checks) {
                if (meas_t_offset == 0) {
                    if (!w->is_x_parity) {
                        sim->m_record.storage[stab_meas_time[w]][0] = 0;
                    }
                } else if (meas_t_offset == n) {
                    // Fix only last round.
                    sim->m_record.storage[n + stab_meas_time[w]][0] = sim->m_record.storage[stab_meas_time[w]][0];
                } else {
                    // Fix last two rounds.
                    sim->m_record.storage[meas_t_offset + stab_meas_time[w]][0] =
                        sim->m_record.storage[meas_t_offset + stab_meas_time[w] - 2*n][0];
                    sim->m_record.storage[meas_t_offset + stab_meas_time[w] - n][0] =
                        sim->m_record.storage[meas_t_offset + stab_meas_time[w] - 2*n][0];
                }
            }
            // Using infected set, identify SWAP LRCs.
            std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*> swap_targets;
            if (perform_swaps) {
                if ((r == circuit_params.rounds - 1) && (flags & EN_STATE_DRAIN)) {
                    // Perform last minute swap LRUs to kill any remaining leakage errors.
                    for (auto pair : swap_set) {
                        swap_targets[pair.first] = pair.second;
                    }
                    // We will simulate "Temporary Stabilizer Extension" on the unlucky data qubit.
                } else {
                    for (auto v : infected) {
                        for (auto w : lattice_graph.adjacency_list(v)) {
                            if (!swap_targets.count(w)) {
                                swap_targets[w] = v;
                                break;
                            }
                        }
                    }
                }
            }
            // Perform operations.
            for (auto v : data_qubits) {
                apply_round_start_error(v->qubit);
            }

            // Update previous syndromes and leakages
            for (uint j = prev_syndrome.size() - 1; j >= 1; j--) {
                prev_syndrome[j] = prev_syndrome[j-1];
                prev_leakages[j] = prev_leakages[j-1];
            }
            prev_syndrome[0] = syndrome;
            prev_leakages[0] = leakages;

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
                syndrome[i] = (sim->m_record.storage[measurement_time][0] ^ sim->m_record.storage[pmt][0])
                                & ~leakages[i];

                measurement_time++;
            }
            if ((flags & EN_TMP_STAB_EXT) && r == circuit_params.rounds-1) {
                // Simulate "Temporary Stabilizer Extension"
                // We apply the same errors as we would at the beginning of a round (depolarizing + leakage)
                if (sim->leakage_table[unlucky_data_qubit->qubit][0]) {
                    // Only situation where the SWAP fails.
                    apply_reset(unlucky_data_qubit->qubit);
                } else {
                    apply_round_start_error(unlucky_data_qubit->qubit);
                }
                // Finally, perform "Drip Testing".
                if (flags & EN_DRIP_TEST) {
                    for (auto v : data_qubits) {
                        auto w = lattice_graph.adjacency_list(v)[0];
                        // It doesn't really matter what happens. We are just simulating the operations.
                        // First, excite v into |2> (or try to).
                        bool orig_x_state = sim->x_table[v->qubit][0];
                        if (sim->x_table[v->qubit][0]) {
                            sim->leakage_table[v->qubit][0] ^= 1;
                        }
                        apply_CX(v->qubit, w->qubit);
                        apply_CX(v->qubit, w->qubit);
                        apply_CX(v->qubit, w->qubit);
                        if (sim->x_table[v->qubit][0]) {
                            sim->leakage_table[v->qubit][0] ^= 1;
                        }
                        // We don't want to measure these qubits directly 
                        // because this will screw up our simulation tables
                        if (sim->leakage_table[v->qubit][0] || (!(rng() & 0x7) && orig_x_state)) {
                            apply_reset(v->qubit);
                        }
                    }
                }
            }

            if (maintain_failure_log) {
                shot_log += "R" + std::to_string(r) + "------------------------\n";
                write_leakage_condition_to_log(shot_log);
                write_aliases_to_log(shot_log, swap_targets);
                write_syndrome_to_log(shot_log, syndrome, leakages);
            }

            meas_t_offset += parity_qubits.size();
        }
        // Perform tail measurements.
        for (auto v : data_qubits) {
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
        syndromes[syndrome_table_index++] |= measure_results[0];
        syndromes[syndrome_table_index++] |= leakage_results[0];
        shots--;

        if (leakage_results[0][base_circuit.count_detectors()] && maintain_failure_log) {
            failure_log += "===============================\n" + shot_log;
        }
    }
    return syndromes;
}

void
Fleece::compute_optimal_swap_set() {
    // Just get a perfect matching between data qubits and {parity qubits U dummy}
    uint number_of_nodes = data_qubits.size() + parity_qubits.size() + 1;
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
        // Now apply leakage error.
        if (postresleak_table.count(qubit)) {
            p = postresleak_table[qubit];
        } else {
            p = circuit_params.get_after_reset_leakage_probability();
            postresleak_table[qubit] = p;
        }

        stim::OperationData lerror_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->LEAKAGE_ERROR(lerror_data);
    }
}

void
Fleece::apply_round_start_error(uint32_t qubit) {
    stim::GateTarget q{qubit};

    fp_t p;
    if (rounddp_table.count(qubit)) {
        p = rounddp_table[qubit];
    } else {
        p = circuit_params.get_before_round_data_depolarization();
        rounddp_table[qubit] = p;
    }
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
Fleece::apply_H(uint32_t qubit) {
    fp_t p = 0;
    stim::GateTarget q{qubit};
    stim::OperationData h_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->H_XZ(h_data);

    if (clifforddp1_table.count(qubit)) {
        p = clifforddp1_table[qubit];
    } else {
        p = circuit_params.get_after_clifford_depolarization(true);
        clifforddp1_table[qubit] = p;
    }
    stim::OperationData dp_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
    sim->DEPOLARIZE1(dp_data);
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
