/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

#include "qontra/ext/qes.h"
#include "qontra/sim/full_system_sim.h"

namespace qontra {
    
inline std::vector<uint32_t>
v32(std::vector<uint64_t> x) {
    return std::vector<uint32_t>(x.begin(), x.end());
}

void
FullSystemSimulator::read_next_instruction(const qes::Program<>& from, program_status_t& st) {
    const qes::Instruction<>& instruction = from.at(st.pc);
    const std::string name = instruction.get_name();

    snapshot();
    // We will have to rollback any instructions if a trial is on a different PC
    // than st.pc. This is just all trials where st.branch_and_wait_reach_map or
    // st.return_if_waiting_trials is set (as these trials are doing nothing).
    stim::simd_bits rollback_pred(st.return_if_waiting_trials);
    for (auto& p : st.branch_and_wait_reach_map) {
        rollback_pred |= p.second;
    }
    
    // Injecting decoherence + dephasing errors.
    if (instruction.has_annotation("timing_error")) {
        inject_timing_error();
    }

    // Identify the type of instruction and act accordingly.
    const isa_data_t& itype = isa_get(instruction);

    if (is_gate(instruction)) {
        elapsed_time += do_gate(instruction);
    } else if (itype == INSTRUCTION_TYPE::EVENT_OR_OBS) {
        create_event_or_obs(instruction);
    } else if (itype == INSTRUCTION_TYPE::MEMORY_SHIFT) {
        int64_t shift = inst.get<int64_t>(0);
        if (instruction.get_name() == "mshift") {
            base_sim->shift_record_by(shift);
            meas_ctr -= shift;
        } else if (instruction.get_name() == "eshift") {
            left_shift(syndrome_table, shift);
        }
    } else if (itype == INSTRUCTION_TYPE::MEMORY_OFFSET) {
        int64_t offset = inst.get<int64_t>(0);
        if (instruction.get_name() == "moffset") {
            meas_offset += offset;
        } else if (instruction.get_name() == "eoffset") {
            event_offset += offset;
        }
    } else if (itype == INSTRUCTION_TYPE::PERMISSIVE_BR_TYPE1) {
        uint64_t target_pc = static_cast<uint64_t>(instruction.get<int64_t>(0));
        std::string r = instruction.get<std::string>(1);

        st.branch_and_wait_reach_map[target_pc] = get_register(r);
    } else if (itype == INSTRUCTION_TYPE::PERMISSIVE_BR_TYPE2) {
        std::string r = instruction.get<std::string>(0);
        st.return_if_waiting_trials |= get_register(r);
    } else if (itype == INSTRUCTION_TYPE::REGISTER_LOGICAL) {
        std::string r1 = instruction.get<std::string>(0);
        if (name == "not") {
            get_register(r1).invert_bits();
        } else {
            std::string r2 = instruction.get<std::string>(1);
            if (name == "or") {
                get_register(r1) |= get_register(r2);
            } else if (name == "and") {
                get_register(r1) &= get_register(r2);
            } else if (name == "xor") {
                get_register(r1) ^= get_register(r2);
            }
        }
    } else if (itype == INSTRUCTION_TYPE::DATA_MOVEMENT) {
        std::string r1 = instruction.get<std::string>(0);
        // We will copy other (well, actually swap) into r1.
        stim::simd_bits<SIMD_WIDTH> other(current_shots);
        if (name == "mov") {
            std::string r2 = instruction.get<std::string>(1);
            other = get_register(r2);
        } else {
            int64_t i = instruction.get<int64_t>(1);
            if (name == "movm") {
                other = base_sim->record_table[i];
            } else if (name == "move") {
                other = syndrome_table[i];
            } else if (name == "movo") {
                other = observable_table[i];
            }
        }
        get_register(r1).swap_with(other);
    }
    // Execute idling errors if necessary.
    if (is_2q_op(instruction)) {
        inject_idling_errors_negative(get_qubits(instruction));
    }
    // Perform rollback.
    rollback_where(rollback_pred);
    // Update PC.
    st.pc++;
    // Update branch data structures. Only update st.branch_and_wait_reach_map,
    // as st.return_if_waiting_trials is microcode-specific.
    st.branch_and_wait_reach_map.erase(st.pc);
}

fp_t
FullSystemSimulator::do_gate(const qes::Instruction<>& instruction, int64_t trial) {
    std::string name = instruction.get_name();
    std::vector<uint64_t> qubits = get_qubits(instruction);
    if (name == "measure") { 
        return do_measurement(instruction, trial);
    } else if (name == "h") {
        sim->H(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("H", v32(qubits));
        }
    } else if (op == "x") {
        sim->X(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("X", v32(qubits));
        }
    } else if (op == "z") {
        sim->Z(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("Z", v32(qubits));
        }
    } else if (op == "cx") {
        sim->CX(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("CX", v32(qubits));
        }
    } else if (op == "reset") {
        sim->R(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("R", v32(qubits));
        }
    } else { return 0.0; }
    // Inject errors.
    std::vector<fp_t> e_array = get_errors(instruction, config.errors);
    if (is_recording_stim_instructions && trial < 0) {
        std::string error_name = apply_x_error_instead_of_depolarizing(instruction)
                                    ? "X_ERROR" : "DEPOLARIZE1";
        sample_circuit.apply_errors(error_name, qubits, e_array, is_2q_op(instruction)); 
    }
    // Now, inject errors. We need to handle the case where trial >= 0
    // differently.
    if (trial >= 0) {
        for (size_t i = 0; i < e_array.size(); i++) {
            fp_t e = e_array[i];
            if (is_2q_op(instruction)) {
                uint64_t q1 = qubits[2*i],
                        q2 = qubits[2*i+1];
                if (sim->get_probability_sample_from_rng() < e) {
                    sim->eDP2(q2, q2, trial);
                }
            } else {
                uint64_t q = qubits[i];
                if (sim->get_probability_sample_from_rng() < e) {
                    if (apply_x_error_instead_of_depolarizing(instruction)) {
                        sim->eX(q, trial);
                    } else {
                        sim->eDP1(q, trial);
                    }
                }
            }
        }
    } else {
        if (is_2q_op(instruction)) {
            base_sim->error_channel<StateSimulator::eDP2>(qubits, e_array);
        } else {
            if (apply_x_error_instead_of_depolarizing(instruction)) {
                base_sim->error_channel<StateSimulator::eX>(qubits, e_array);
            } else {
                base_sim->error_channel<StateSimulator::eDP1>(qubits, e_array);
            }
        }
    }
    // Return instruction latency.
    return get_max_latency(instruction, config.timing);
}

fp_t
FullSystemSimulator::do_measurement(const qes::Instruction<>& instruction, int64_t trial) {
    std::vector<uint64_t> qubits = get_qubits(instruction);

    std::vector<fp_t> m1w0_array, m0w1_array;
    for (size_t i = 0; i < qubits.size(); i++) {
        uint64_t q = qubits[i];
        m1w0_array.push_back(config.errors.m1w0[q]);
        m0w1_array.push_back(config.errors.m0w1[q]);
        // Stim doesn't have a way to implement biased measurements.
        // So, we just take the mean.
        if (is_recording_stim_instructions && trial < 0) {
            fp_t emean = 0.5 * (config.errors.m1w0[q]+config.errors.m0w1[q]);
            sample_circuit.safe_append_ua("X_ERROR", {q}, emean);
        }
    }
    // Perform the measurement.
    sim->M(operands, m1w0_array, m0w1_array, meas_ctr+meas_offset, trial);
    if (is_recording_stim_instructions && trial < 0) {
        sample_circuit.safe_append_u("M", operands);
    }
    // We will only update the counter if this is a common measurement.
    if (trial < 0)  meas_ctr += operands.size();
    return get_max_latency(instruction, config.timing);
}


}   // qontra
