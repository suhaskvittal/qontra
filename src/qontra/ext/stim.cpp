/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include "qontra/ext/stim.h"
#include "qontra/ext/qes.h"

#include <algorithm>

#include <math.h>

static const uint32_t STIM_INSTRUCTION_LEN_LIMIT = 64;

namespace qontra {

inline std::vector<uint32_t>
v32(std::vector<uint64_t> arr) {
    return std::vector<uint32_t>(arr.begin(), arr.end());
}

DetailedStimCircuit::DetailedStimCircuit()
    :stim::Circuit(),
    number_of_colors_in_circuit(0),
    detector_base_map(),
    detector_color_map(),
    flag_detectors()
{}

DetailedStimCircuit::DetailedStimCircuit(const stim::Circuit& other)
    :stim::Circuit(other),
    number_of_colors_in_circuit(0),
    detector_base_map(),
    detector_color_map(),
    flag_detectors()
{}

DetailedStimCircuit::DetailedStimCircuit(const DetailedStimCircuit& other)
    :stim::Circuit(other),
    number_of_colors_in_circuit(other.number_of_colors_in_circuit),
    detector_base_map(other.detector_base_map),
    detector_color_map(other.detector_color_map),
    flag_detectors(other.flag_detectors)
{}

DetailedStimCircuit
DetailedStimCircuit::from_qes(
        const qes::Program<>& program,
        const ErrorTable& errors,
        const TimeTable& timing,
        fp_t fix_timing_error_as_depolarizing_error) 
{
    DetailedStimCircuit circuit;

    std::set<int> all_colors;

    std::vector<uint64_t> all_qubits;
    get_number_of_qubits(program, all_qubits);

    fp_t elapsed_time = 0.0;
    size_t meas_ctr = 0;
    for (const qes::Instruction<>& inst : program) {
        std::string name = inst.get_name();
        std::vector<uint64_t> qubits = get_qubits(inst);
        // Read annotations of instruction.
        bool has_timing_error = inst.has_annotation("timing_error");
        bool is_error_free = inst.has_annotation("no_error");
        bool takes_no_time = inst.has_annotation("no_tick") || is_instantaneous(inst);
        // Check if we need to inject a timing error.
        if (has_timing_error) {
            if (fix_timing_error_as_depolarizing_error >= 0.0) {
                circuit.safe_append_ua("DEPOLARIZE1", v32(all_qubits), fix_timing_error_as_depolarizing_error);
            } else {
                for (size_t i = 0; i < qubits.size(); i++) {
                    fp_t t1 = timing.t1.at(qubits[i]),
                         t2 = timing.t2.at(qubits[i]);
                    fp_t e_ad = 0.25 * (1 - exp(-elapsed_time/t1));
                    fp_t e_pd = 0.5 * (1 - exp(-elapsed_time/t2));

                    uint32_t _q = static_cast<uint32_t>(qubits[i]);

                    circuit.safe_append_ua("X_ERROR", {_q}, e_ad);
                    circuit.safe_append_ua("Y_ERROR", {_q}, e_ad);
                    circuit.safe_append_ua("Z_ERROR", {_q}, e_pd-e_ad);
                }
            }
            elapsed_time = 0.0;
        }
        if (isa_get(inst).error_precedes_op() && !is_error_free) {
            std::vector<fp_t> e_array = get_errors(inst, errors);
            std::string error_name = isa_get(inst).apply_x_error_instead_of_depolarizing() 
                                        ? "X_ERROR" : "DEPOLARIZE1";
            circuit.apply_errors(error_name, qubits, e_array, is_2q_gate(inst));
        }
        // Convert the qes::Instruction to stim if there's an equivalent.
        if (name == "h") {
            circuit.safe_append_u("H", v32(qubits));
        } else if (name == "x") {
            circuit.safe_append_u("X", v32(qubits));
        } else if (name == "z") {
            circuit.safe_append_u("Z", v32(qubits));
        } else if (name == "cx") {
            circuit.safe_append_u("CX", v32(qubits));
        } else if (name == "cz") {
            circuit.safe_append_u("CZ", v32(qubits));
        } else if (name == "measure") {
            circuit.safe_append_u("M", v32(qubits));
            meas_ctr += qubits.size();
        } else if (name == "reset") {
            circuit.safe_append_u("R", v32(qubits));
        } else if (name == "event") {
            std::vector<uint32_t> offsets;
            int64_t detection_event = inst.get<int64_t>(0);
            // Note: first operand is not useful to Stim.
            for (size_t i = 1; i < inst.get_number_of_operands(); i++) {
                uint32_t off = static_cast<uint32_t>(meas_ctr - inst.get<int64_t>(i));
                offsets.push_back(stim::TARGET_RECORD_BIT | off);
            }
            // Check annotations and property map for any additional data.
            int color_id = 0;
            if (inst.has_property("base")) {
                uint64_t base = static_cast<uint64_t>(inst.get_property<int64_t>("base"));
                circuit.detector_base_map[detection_event] = base;
            }
            if (inst.has_property("color")) {
                color_id = static_cast<int>(inst.get_property<int64_t>("color"));
                all_colors.insert(color_id);
                circuit.detector_color_map[detection_event] = color_id;
            }
            if (inst.has_annotation("flag")) {
                circuit.flag_detectors.insert(detection_event);
            }
            // Append instruction
            const std::vector<double> coord{
                static_cast<double>(detection_event),
                0.0,
                0.0,
                static_cast<double>(color_id)
            };
            circuit.safe_append_u("DETECTOR", offsets, coord);
        } else if (name == "obs") {
            std::vector<uint32_t> offsets;
            int64_t obs = inst.get<int64_t>(0);
            for (size_t i = 1; i < inst.get_number_of_operands(); i++) {
                uint32_t off = static_cast<uint32_t>(meas_ctr - inst.get<int64_t>(i));
                offsets.push_back(stim::TARGET_RECORD_BIT | off);
            }
            circuit.safe_append_ua("OBSERVABLE_INCLUDE", offsets, static_cast<double>(obs));
        } else if (name == "mshift") {
            meas_ctr -= inst.get<int64_t>(0);
        }
        // Apply any errors.
        // First do gate errors.
        if (!is_error_free && !isa_get(inst).error_precedes_op()) {
            std::vector<fp_t> e_array = get_errors(inst, errors);
            std::string error_name = isa_get(inst).apply_x_error_instead_of_depolarizing() 
                                        ? "X_ERROR" : "DEPOLARIZE1";
            circuit.apply_errors(error_name, qubits, e_array, is_2q_gate(inst));
        }
        // Now do idling errors and update elapsed_time.
        if (!takes_no_time) {
            // Use uint32_t as Stim uses that (saves a cast).
            if (is_2q_gate(inst) && !is_error_free) {
                for (uint64_t q : all_qubits) {
                    if (std::find(qubits.begin(), qubits.end(), q) == qubits.end()) {
                        // This will have an idling error.
                        fp_t e = errors.idling.at(q);
                        circuit.safe_append_ua("DEPOLARIZE1", {static_cast<uint32_t>(q)}, e);
                    }
                }
            }
            elapsed_time += get_max_latency(inst, timing);
        }
    }
    circuit.number_of_colors_in_circuit = all_colors.size();
    return circuit;
}

}   // qontra
