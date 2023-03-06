/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#include "fleece.h"

namespace qrc {

//#define FLEECE_DEBUG

template <typename T>
using PTR = stim::ConstPointerRange<T>;

Fleece::Fleece(const stim::CircuitGenParameters& params,
                std::mt19937_64& rng,
                char reset_basis,
                char output_basis)
    :circuit_params(params),
    reset_basis(reset_basis),
    output_basis(output_basis),
    base_circuit(stim::generate_surface_code_circuit(params).circuit),
    sim(nullptr),
    lattice_graph(),
    data_qubits(),
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
}

Fleece::~Fleece() {
    delete sim;
}

stim::simd_bit_table
Fleece::create_syndromes(uint64_t shots, uint disable_leakage_at_round) {
    const size_t num_results = base_circuit.count_detectors() + base_circuit.count_observables() + 1;
    stim::simd_bit_table syndromes(2*shots, num_results);
    syndromes.clear();

    uint64_t syndrome_table_index = 0;
    while (shots) {
        sim->reset_all();
        sim->leakage_enabled = true;
        //
        // Do first round -- special.
        //
        // Reset all qubits.
        uint32_t measurement_time = 0;
        std::vector<uint8_t> syndrome(parity_qubits.size(), 0);
        std::vector<uint8_t> leakages(parity_qubits.size(), 0);
        std::vector<fleece::LatticeGraph::Vertex*> acting_parity_qubits(parity_qubits.size());

        std::map<fleece::LatticeGraph::Vertex*, uint32_t> stab_meas_time;
            
        for (auto v : data_qubits) {
            apply_reset(v->qubit, false);
            apply_round_start_error(v->qubit);
        }
        for (auto v : parity_qubits) {
            apply_reset(v->qubit, false);
        }
        // Perform stabilizer measurements.
        for (auto v : parity_qubits) {
            if (v->is_x_parity) {
                apply_H(v->qubit);
                for (auto w : lattice_graph.adjacency_list(v)) {
                    apply_CX(v->qubit, w->qubit);
                }
                apply_H(v->qubit);
            } else {
                for (auto w : lattice_graph.adjacency_list(v)) {
                    apply_CX(w->qubit, v->qubit);
                }
            }
            apply_measure(v->qubit);
            apply_reset(v->qubit);
            if (v->is_x_parity) {
                syndrome[measurement_time] = 0;
            } else {
                // Only record as error if parity qubit isn't leaked.
                syndrome[measurement_time] = sim->m_record.storage[measurement_time][0]
                                                & ~sim->leak_record.storage[measurement_time][0];
            }
            leakages[measurement_time] = sim->leak_record.storage[measurement_time][0];
            acting_parity_qubits[measurement_time] = v;
            stab_meas_time[v] = measurement_time;
            measurement_time++;
        }
        // Now, we need to simulate the circuit operation by operation, round by round.
        uint32_t meas_t_offset = 0;
        for (uint r = 1; r < circuit_params.rounds; r++) {
            if (r == disable_leakage_at_round) {
                sim->toggle_leakage();
            }
            std::set<fleece::LatticeGraph::Vertex*> infected;   // Set of vertices with potential leakages.
            // Check for problems in the prior round.
            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                if (leakages[i]) {
                    auto v = acting_parity_qubits[i];
                    std::cout << "leak detected on qubit " << v->qubit << "\n";
                    // If v is a data qubit, then only the adjacent stabilizers were affected.
                    // So, clear the measurements.
                    //
                    // If v is a parity qubit, then we must also match flipped stabilizers to
                    // this leaked stabilizer measurement.
                    //
                    // TODO: account for measurement/CNOT errors on flipped stabilizers.
                    if (v->is_data) {
                        for (auto w : lattice_graph.adjacency_list(v)) {
                            if (meas_t_offset == 0) {
                                if (!w->is_x_parity) {
                                    sim->m_record.storage[stab_meas_time[w]][0] = 0;
                                }
                            } else {
                                sim->m_record.storage[meas_t_offset + stab_meas_time[w]][0] =
                                    sim->m_record.storage[meas_t_offset + stab_meas_time[w] - parity_qubits.size()][0];
                            }
                        }
                    } else {
                        for (auto w : parity_qubits) {
                            uint32_t mt = stab_meas_time[w];
                            if (!syndrome[mt]) {
                                continue;
                            }
                            auto common = lattice_graph.get_common_neighbors(v, w);
                            for (auto u : common) {
                                if (w->is_x_parity) {
                                    sim->z_table[u->qubit][0] ^= 1;
                                } else {
                                    sim->x_table[u->qubit][0] ^= 1;
                                }
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
            // Using infected set, identify SWAP LRCs.
            std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*> swap_targets;
            for (auto v : infected) {
                for (auto w : lattice_graph.adjacency_list(v)) {
                    if (!swap_targets.count(w)) {
                        swap_targets[w] = v;
                        break;
                    }
                }
            }
            // Perform operations.
            for (auto v : data_qubits) {
                apply_round_start_error(v->qubit);
            }

            for (uint32_t i = 0; i < parity_qubits.size(); i++) {
                auto v = parity_qubits[i];
                if (v->is_x_parity) {
                    apply_H(v->qubit);
                    for (auto w : lattice_graph.adjacency_list(v)) {
                        apply_CX(v->qubit, w->qubit);
                    }
                    apply_H(v->qubit);
                } else {
                    for (auto w : lattice_graph.adjacency_list(v)) {
                        apply_CX(w->qubit, v->qubit);
                    }
                }
                
                if (swap_targets.count(v)) {
                    auto _v = swap_targets[v];
                    apply_SWAP(v->qubit, _v->qubit);
                    apply_measure(_v->qubit);
                    apply_reset(_v->qubit);
                    // Don't apply a depolarizing error if this is the last syndrome extraction round.
                    apply_SWAP(v->qubit, _v->qubit, (r != circuit_params.rounds - 1));

                    acting_parity_qubits[i] = _v;
                } else {
                    apply_measure(v->qubit);
                    apply_reset(v->qubit);

                    acting_parity_qubits[i] = v;
                }
                uint32_t pmt = measurement_time - parity_qubits.size();
                syndrome[i] = (sim->m_record.storage[measurement_time][0] ^ sim->m_record.storage[pmt][0])
                                & ~sim->leak_record.storage[measurement_time][0];
                leakages[i] = sim->leak_record.storage[measurement_time][0];

                measurement_time++;
            }
            meas_t_offset += parity_qubits.size();
        }
        // Perform tail measurements.
        for (auto v : data_qubits) {
            apply_measure(v->qubit, false);
        }
        // Everything is done, extract the syndrome.
        size_t num_results = base_circuit.count_detectors() + base_circuit.count_observables() + 1;
        stim::simd_bit_table measure_results(num_results, 1);
        stim::simd_bit_table leakage_results(num_results, 1);
        auto det_obs = stim::DetectorsAndObservables(base_circuit);
        stim::read_from_sim(*sim, det_obs, false, true, true,
                                measure_results, leakage_results);
        measure_results = measure_results.transposed();
        leakage_results = leakage_results.transposed();
#ifdef FLEECE_DEBUG
        auto msyn = _to_vector(measure_results[0], base_circuit.count_detectors(), base_circuit.count_observables());
        auto lsyn = _to_vector(leakage_results[0], base_circuit.count_detectors(), base_circuit.count_observables());

        const uint d = circuit_params.distance;
        const uint n_detectors = (d*d - 1) >> 1;
        std::cout << "================================\n";
        uint r = 0;
        for (uint i = 0; i < base_circuit.count_detectors(); i += n_detectors) {
            std::cout << r << "\t";
            for (uint j = 0; j < n_detectors; j++) {
                if (lsyn[i+j]) {
                    std::cout << "L";
                } else if (msyn[i+j]) {
                    std::cout << "!";
                } else {
                    std::cout << ".";
                }
            }
            std::cout << "\n";
            r++;
        }
        const uint last_index = base_circuit.count_detectors();
        if (lsyn[last_index]) {
            std::cout << "L\tL\n";
        } else if (msyn[last_index]) {
            std::cout << "L\t!\n";
        } else {
            std::cout << "L\t.\n";
        }
#endif
        syndromes[syndrome_table_index++] |= measure_results[0];
        syndromes[syndrome_table_index++] |= leakage_results[0];
        shots--;
    }
    return syndromes;
}

void
Fleece::apply_reset(int32_t qubit, bool add_error) {
    stim::GateTarget q{qubit};
    stim::OperationData reset_data{PTR<double>(), PTR<stim::GateTarget>(&q)};
    sim->reset_z(reset_data);

    if (add_error) {
        fp_t p;
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
Fleece::apply_round_start_error(int32_t qubit) {
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
Fleece::apply_H(int32_t qubit) {
    stim::GateTarget q{qubit};
    stim::OperationData h_data{PTR<double>(), PTR<stim::GateTarget>(&q)};
    sim->H_XZ(h_data);

    fp_t p;
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
Fleece::apply_CX(int32_t qubit1, int32_t qubit2) {
    stim::GateTarget q1{qubit1};
    stim::GateTarget q2{qubit2};
    std::vector<stim::GateTarget> gate_targets{q1, q2};
    stim::OperationData cx_data{PTR<double>(), PTR<stim::GateTarget>(gate_targets)};
    sim->ZCX(cx_data);

    fp_t p;
    std::vector<int32_t> q1_q2{qubit1, qubit2};
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

void
Fleece::apply_measure(int32_t qubit, bool add_error) {
    stim::GateTarget q{qubit};
    if (add_error) {
        fp_t p;
        if (premeasflip_table.count(qubit)) {
            p = premeasflip_table[qubit];
        } else {
            p = circuit_params.get_before_measure_flip_probability();
            premeasflip_table[qubit] = p;
        }
        stim::OperationData flip_data{PTR<double>(&p), PTR<stim::GateTarget>(&q)};
        sim->X_ERROR(flip_data);
    }

    stim::OperationData meas_data{PTR<double>(), PTR<stim::GateTarget>(&q)};
    sim->measure_z(meas_data);
}

void
Fleece::apply_SWAP(int32_t qubit1, int32_t qubit2, bool add_error) {
    stim::GateTarget q1{qubit1};
    stim::GateTarget q2{qubit2};
    std::vector<stim::GateTarget> gate_targets{q1, q2};
    stim::OperationData swap_data{PTR<double>(), PTR<stim::GateTarget>(gate_targets)};
    sim->SWAP(swap_data);

    if (add_error) {
        fp_t p;
        std::vector<int32_t> q1_q2{qubit1, qubit2};
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
