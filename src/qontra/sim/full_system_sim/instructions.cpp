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
    INSTRUCTION_TYPE itype = isa_get(instruction).inst_type;

    if (is_gate(instruction)) {
        fp_t lat = do_gate(instruction);
        // Update timing only where rollback_pred = 0.
        elapsed_time += lat;
        for (uint64_t t = 0; t < current_shots; t++) {
            if (rollback_pred[t]) shot_time_delta_map[t] -= lat;
        }
    } else if (itype == INSTRUCTION_TYPE::EVENT_OR_OBS) {
        create_event_or_obs(instruction);
    } else if (itype == INSTRUCTION_TYPE::MEMORY_SHIFT) {
        int64_t shift = instruction.get<int64_t>(0);
        if (instruction.get_name() == "mshift") {
            base_sim->shift_record_by(shift);
            meas_ctr -= shift;
        } else if (instruction.get_name() == "eshift") {
            left_shift(syndrome_table, shift);
        }
    } else if (itype == INSTRUCTION_TYPE::MEMORY_OFFSET) {
        int64_t offset = instruction.get<int64_t>(0);
        if (instruction.get_name() == "moffset") {
            meas_offset += offset;
        } else if (instruction.get_name() == "eoffset") {
            event_offset += offset;
        }
    } else if (itype == INSTRUCTION_TYPE::PERMISSIVE_BR_TYPE1) {
        uint64_t target_pc = static_cast<uint64_t>(instruction.get<int64_t>(0));
        std::string r = instruction.get<std::string>(1);

        st.branch_and_wait_reach_map.emplace(target_pc, get_register(r));
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
    if (is_2q_gate(instruction)) {
        inject_idling_error_negative(get_qubits(instruction));
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
        base_sim->H(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("H", v32(qubits));
        }
    } else if (name == "x") {
        base_sim->X(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("X", v32(qubits));
        }
    } else if (name == "z") {
        base_sim->Z(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("Z", v32(qubits));
        }
    } else if (name == "s") {
        base_sim->S(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("S", v32(qubits));
        }
    } else if (name == "sdg") {
        base_sim->S(qubits, trial);
        base_sim->Z(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("S_DAG", v32(qubits));
        }
    } else if (name == "cx") {
        base_sim->CX(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("CX", v32(qubits));
        }
    } else if (name == "reset") {
        base_sim->R(qubits, trial);
        if (is_recording_stim_instructions && trial < 0) {
            sample_circuit.safe_append_u("R", v32(qubits));
        }
    } else { return 0.0; }
    // Inject errors.
    std::vector<fp_t> e_array = get_errors(instruction, config.errors);
    if (is_recording_stim_instructions && trial < 0) {
        std::string error_name = isa_get(instruction).apply_x_error_instead_of_depolarizing() 
                                    ? "X_ERROR" : "DEPOLARIZE1";
        sample_circuit.apply_errors(error_name, qubits, e_array, is_2q_gate(instruction)); 
    }
    // Now, inject errors. We need to handle the case where trial >= 0
    // differently.
    if (trial >= 0) {
        for (size_t i = 0; i < e_array.size(); i++) {
            fp_t e = e_array[i];
            if (is_2q_gate(instruction)) {
                uint64_t q1 = qubits[2*i],
                        q2 = qubits[2*i+1];
                if (base_sim->get_probability_sample_from_rng() < e) {
                    base_sim->eDP2(q1, q2, trial);
                }
            } else {
                uint64_t q = qubits[i];
                if (base_sim->get_probability_sample_from_rng() < e) {
                    if (isa_get(instruction).apply_x_error_instead_of_depolarizing()) {
                        base_sim->eX(q, trial);
                    } else {
                        base_sim->eDP1(q, trial);
                    }
                }
            }
        }
    } else {
        if (is_2q_gate(instruction)) {
            base_sim->error_channel<&StateSimulator::eDP2>(qubits, e_array);
        } else {
            if (isa_get(instruction).apply_x_error_instead_of_depolarizing()) {
                base_sim->error_channel<&StateSimulator::eX>(qubits, e_array);
            } else {
                base_sim->error_channel<&StateSimulator::eDP1>(qubits, e_array);
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
            uint32_t _q = static_cast<uint32_t>(q);
            fp_t emean = 0.5 * (config.errors.m1w0[q]+config.errors.m0w1[q]);
            sample_circuit.safe_append_ua("X_ERROR", {_q}, emean);
        }
    }
    // Perform the measurement.
    base_sim->M(qubits, m1w0_array, m0w1_array, meas_ctr+meas_offset, trial);
    if (is_recording_stim_instructions && trial < 0) {
        sample_circuit.safe_append_u("M", v32(qubits));
    }
    // We will only update the counter if this is a common measurement.
    if (trial < 0)  meas_ctr += qubits.size();
    return get_max_latency(instruction, config.timing);
}

void
FullSystemSimulator::create_event_or_obs(const qes::Instruction<>& instruction) {
    std::string name = instruction.get_name();
    int64_t index = instruction.get<int64_t>(0);
    if (name == "event") index += event_offset;
    
    stim::simd_bits_range_ref<SIMD_WIDTH> w = 
        name == "event" ? syndrome_table[index] : observable_table[index];
    w.clear();
    uint64_t& max_index = name == "event" ? max_event_written_to : max_obs_written_to;
    for (size_t i = 1; i < instruction.get_number_of_operands(); i++) {
        int64_t k = instruction.get<int64_t>(i);
        w ^= base_sim->record_table[k];
    }
    max_index = std::max(static_cast<uint64_t>(index)+1, max_index);

    if (is_recording_stim_instructions) {
        std::vector<uint32_t> offsets;
        for (size_t i = 1; i < instruction.get_number_of_operands(); i++) {
            int64_t k = instruction.get<int64_t>(i);
            offsets.push_back(stim::TARGET_RECORD_BIT | static_cast<uint32_t>(meas_ctr+meas_offset-k));
        }
        const std::vector<fp_t> coord{static_cast<fp_t>(index)};

        std::string stim_inst_name = name == "event" ? "DETECTOR" : "OBSERVABLE_INCLUDE";
        sample_circuit.safe_append_u(stim_inst_name, offsets, coord);
    }
}


}   // qontra
