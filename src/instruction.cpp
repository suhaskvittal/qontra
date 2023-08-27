/*
 *  author: Suhas Vittal date:   24 June 2023
 * */

#include "instruction.h"

namespace qontra {

std::vector<uint>
Instruction::get_qubit_operands() const {
    // Note: this is for qubit operands that are physically interacted with!
    //
    // Do not include instructions like lockq or unlockq which do not interact
    // with qubits directly.
    if (ONLY_HAS_QUBIT_OPERANDS.count(name)) {
        return operands;
    } else if (name == "mrc") {
        return std::vector<uint>(operands.begin()+1, operands.end());
    } else {
        return std::vector<uint>();
    }
}

uint
get_number_of_qubits(const schedule_t& sch) {
    uint max_q = 0;
    for (const auto& inst : sch) {
        for (uint x : inst.operands) {
            max_q = x > max_q ? x : max_q;
        }
    }
    return max_q+1;
}

std::string
schedule_to_text(const schedule_t& sch) {
    std::string out;
    for (const auto& inst : sch) {
        out += inst.str() + ";\n";
    }
    return out;
}

stim::Circuit
schedule_to_stim(const schedule_t& sch, ErrorTable& errors, TimeTable& timing) {
    uint n = get_number_of_qubits(sch);

    stim::Circuit circuit;
    fp_t time = 0.0;
    // The instructions below can only be Pauli operators, CX,
    // or the following:
    //   measure, reset, event, obs, nop
    uint n_meas = 0;
    for (const auto& inst : sch) {
        std::vector<uint> operands = inst.operands;

        bool inject_op_error = !inst.annotations.count(Annotation::no_error);
        bool inject_timing_error = inst.annotations.count(Annotation::inject_timing_error);
        bool operation_takes_time = !inst.annotations.count(Annotation::no_tick);
        bool is_2q_op = IS_2Q_OPERATOR.count(inst.name);
        if (inst.name == "measure") {
            // Add X error before.
            if (inject_op_error) {
                for (uint x : operands) {
                    fp_t e = 0.5*(errors.m1w0[x] + errors.m0w1[x]);
                    if (e > 0)  circuit.append_op("X_ERROR", {x}, e);
                }
            }
            circuit.append_op("M", operands);
            n_meas += operands.size();
        } else if (inst.name == "h") {
            circuit.append_op("H", operands);
        } else if (inst.name == "x") {
            circuit.append_op("X", operands);
        } else if (inst.name == "z") {
            circuit.append_op("Z", operands);
        } else if (inst.name == "cx") {
            circuit.append_op("CX", operands);
        } else if (inst.name == "nop") {
            circuit.append_op("TICK", {});
        } else if (inst.name == "reset") {
            circuit.append_op("R", operands);
        } else if (inst.name == "event") {
            std::vector<uint> offsets;
            for (uint i = 1; i < operands.size(); i++) {
                offsets.push_back(stim::TARGET_RECORD_BIT | (n_meas - operands[i]));
            }
            circuit.append_op("DETECTOR", offsets);
            continue;
        } else if (inst.name == "obs") {
            std::vector<uint> offsets;
            for (uint i = 1; i < operands.size(); i++) {
                offsets.push_back(stim::TARGET_RECORD_BIT | (n_meas - operands[i]));
            }
            circuit.append_op("OBSERVABLE_INCLUDE", offsets, operands[0]);
            continue;
        } else { continue; }    // Ignore all other instructions.
        
        // Inject timing errors if requested.
        if (inject_timing_error) {
            for (uint x = 0; x < n; x++) {
                fp_t t1 = timing.t1[x];
                fp_t t2 = timing.t2[x];

                fp_t e_ad = 0.25*(1 - exp(-time/t1));
                fp_t e_pd = 0.5*(1 - exp(-time/t2));

                circuit.append_op("X_ERROR", {x}, e_ad);
                circuit.append_op("Y_ERROR", {x}, e_ad);
                circuit.append_op("Z_ERROR", {x}, e_pd-e_ad);
            }
            time = 0;
        }
        // Update operation latency.
        if (operation_takes_time) {
            fp_t max_t = 0.0;
            if (is_2q_op) {
                for (uint i = 0; i < operands.size(); i += 2) {
                    fp_t x = operands[i], y = operands[i+1];
                    fp_t t = timing.op2q[inst.name][std::make_pair(x, y)];
                    max_t = t > max_t ? t : max_t;
                }
            } else {
                for (uint x : operands) {
                    fp_t t = timing.op1q[inst.name][x];
                    max_t = t > max_t ? t : max_t;
                }
            }
            time += max_t;
        }
        // Add operation error.
        if (inst.name != "measure" && inject_op_error) {
            if (is_2q_op) {
                for (uint i = 0; i < operands.size(); i += 2) {
                    uint x = operands[i], y = operands[i+1];
                    auto x_y = std::make_pair(x, y);
                    fp_t dp = errors.op2q[inst.name][x_y];
                    fp_t lt = errors.op2q_leakage_transport[inst.name][x_y];
                    fp_t li = errors.op2q_leakage_injection[inst.name][x_y];

                    if (lt > 0) circuit.append_op("L_TRANSPORT", {x, y}, lt);
                    if (li > 0) circuit.append_op("L_ERROR", {x, y}, lt);
                    if (dp > 0) circuit.append_op("DEPOLARIZE2", {x, y}, dp);
                }
            } else {
                for (uint x : operands) {
                    fp_t e = errors.op1q[inst.name][x];
                    if (e > 0) {
                        if (inst.name == "reset") {
                            circuit.append_op("X_ERROR", {x}, e);
                        } else {
                            circuit.append_op("DEPOLARIZE1", {x}, e);
                        }
                    }
                }
            }
        }
    }
    return circuit;
}

}   // qontra
