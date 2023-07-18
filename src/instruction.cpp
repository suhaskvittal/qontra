/*
 *  author: Suhas Vittal
 *  date:   24 June 2023
 * */

#include "instruction.h"

namespace qontra {

std::vector<uint>
Instruction::get_qubit_operands() const {
    if (ONLY_HAS_QUBIT_OPERANDS.count(name)) {
        return operands;
    } else if (name == "mrc") {
        return std::vector<uint>(operands.begin()+1, operands.end());
    } else {
        return std::vector<uint>();
    }
}

std::string
schedule_to_text(const schedule_t& sch) {
    std::string out;
    for (const auto& inst : sch) {
        out += inst.str() + ";\n";
    }
    return out;
}

schedule_t
schedule_from_stim(const stim::Circuit& circ) {
    schedule_t sch;

    uint64_t mctr = 0;  // Measurement, event, and observable counters.
    uint64_t ectr = 0;
    uint64_t octr = 0;
    auto cb = [&] (const stim::Operation& op)
    {
        auto gate = op.gate;
        std::string opname(gate->name);
        std::vector<uint> operands;
        // Parse operands:
        //  We need to handle detector and observable instructions
        //  differently, as they are a record lookback.
        //
        //  We also need to handle measurement instructions as well.
        if (opname == "DETECTOR")               operands.push_back(ectr++);
        if (opname == "OBSERVABLE_INCLUDE")     operands.push_back(octr++);
        if (opname == "M" || opname == "MR")    operands.push_back(mctr);
        for (auto target : op.target_data.targets) {
            uint32_t x = target.data;
            uint32_t v = x & stim::TARGET_VALUE_MASK;
            if (opname == "DETECTOR" || opname == "OBSERVABLE_INCLUDE") {
                operands.push_back(mctr - v);
            } else {
                operands.push_back(v);
            }
        }
        if (opname == "M" || opname == "MR")    mctr += operands.size()-1;
        // Convert stim instruction to ISA instruction
        std::string name;
        if (opname == "X")                  name = "x";
        else if (opname == "Z")             name = "z";
        else if (opname == "CX")            name = "cx";
        else if (opname == "H")             name = "h";
        else if (opname == "R")             name = "reset";
        else if (opname == "S")             name = "s";
        else if (opname == "M")             name = "mrc";
        else if (opname == "DETECTOR")      name = "event";
        else if (opname == "OBSERVABLE_INCLUDE") {
            name = "obs";
        } else if (opname == "MR") {
            sch.push_back({"mrc", operands, {}});
            operands.erase(operands.begin());
            name = "reset";
        } else                              name = "nop";
        if (name == "nop")  return;
        
        sch.push_back({name, operands, {}});
    };
    circ.for_each_operation(cb);
    sch = relabel_operands(sch);
    return sch;
}

schedule_t
relabel_operands(const schedule_t& sch) {
    schedule_t new_sch;

    std::map<uint, uint> operand_map;
    uint k = 0;
    for (const auto& inst : sch) {
        Instruction new_inst;
        new_inst.name = inst.name;
        if (ONLY_HAS_QUBIT_OPERANDS.count(inst.name)) {
            for (uint i : inst.operands) {
                if (!operand_map.count(i))  operand_map[i] = k++;
                uint x = operand_map[i];
                new_inst.operands.push_back(x);
            }
        } else if (inst.name == "mrc") {
            new_inst.operands.push_back(inst.operands[0]);
            for (uint j = 1; j < inst.operands.size(); j++) {
                uint i = inst.operands[j];
                if (!operand_map.count(i))  operand_map[i] = k++;
                uint x = operand_map[i];
                new_inst.operands.push_back(x);
            }
        } else {
            new_inst.operands = inst.operands;
        }
        new_sch.push_back(new_inst);
    }

    return new_sch;
}

stim::Circuit
fast_convert_to_stim(const schedule_t& prog, ErrorTable& errors, TimeTable& timing) {
    stim::Circuit circuit;
    uint mctr = 0;
    for (const auto& inst : prog) {
        // If it is a measurement, inject a measurement error here.
        if (inst.name == "mrc") {
            for (uint i : inst.get_qubit_operands()) {
                circuit.append_op("X_ERROR", {i}, errors.op1q["m"][i]);
            }
        }

        if (inst.name == "x") {
            circuit.append_op("X", inst.get_qubit_operands());
        } else if (inst.name == "z") {
            circuit.append_op("Z", inst.get_qubit_operands());
        } else if (inst.name == "s") {
            circuit.append_op("S", inst.get_qubit_operands());
        } else if (inst.name == "cx") {
            circuit.append_op("CX", inst.get_qubit_operands());
        } else if (inst.name == "mrc") {
            circuit.append_op("M", inst.get_qubit_operands());
            mctr += inst.get_qubit_operands().size();
        } else if (inst.name == "reset") {
            circuit.append_op("R", inst.get_qubit_operands());
        } else if (inst.name == "event") {
            std::vector<uint32_t> components;
            for (uint i = 1; i < inst.operands.size(); i++) {
                components.push_back(stim::TARGET_RECORD_BIT | (mctr - inst.operands[i]));
            }
            circuit.append_op("DETECTOR", components);
        } else if (inst.name == "obs") {
            std::vector<uint32_t> components;
            for (uint i = 1; i < inst.operands.size(); i++) {
                components.push_back(stim::TARGET_RECORD_BIT | (mctr - inst.operands[i]));
            }
            circuit.append_op("OBSERVABLE_INCLUDE", components, inst.operands[0]);
        }
        
        // Add other operation errors and decoherence/dephasing errors
        std::map<uint, fp_t> qubit_to_delay;
        if (inst.name == "cx") {
            for (uint i = 0; i < inst.operands.size(); i += 2) {
                uint j1 = inst.operands[i];
                uint j2 = inst.operands[i+1];
                auto j1_j2 = std::make_pair(j1, j2);

                fp_t dp2 = errors.op2q["cx"][j1_j2];
                fp_t li = errors.op2q_leakage_injection["cx"][j1_j2];
                fp_t lt = errors.op2q_leakage_transport["cx"][j1_j2];

                circuit.append_op("L_TRANSPORT", {j1, j2}, lt);
                circuit.append_op("L_ERROR", {j1, j2}, li);
                circuit.append_op("DEPOLARIZE2", {j1, j2}, dp2);

                fp_t cx_t = timing.op2q["cx"][j1_j2];
                qubit_to_delay[j1] = cx_t;
                qubit_to_delay[j2] = cx_t;
            }
        } else {
            if (inst.name != "mrc") {
                if (!errors.op1q.count(inst.name))   continue;
                for (uint i : inst.operands) {
                    fp_t e = errors.op1q[inst.name][i];
                    if (inst.name == "reset") {
                        circuit.append_op("X_ERROR", {i}, e);
                    } else {
                        circuit.append_op("DEPOLARIZE1", {i}, e);
                    }
                    fp_t t = timing.op1q[inst.name][i];
                    qubit_to_delay[i] = t;
                }
            } else {
                for (uint i : inst.get_qubit_operands()) {
                    fp_t t = timing.op1q["m"][i];
                    qubit_to_delay[i] = t;
                }
            }
        }

        for (auto pair : qubit_to_delay) {
            uint i = pair.first;
            fp_t t = pair.second;
            fp_t mean_t = (timing.t1[i] + timing.t2[i])*0.5;
            fp_t e = 1 - exp(-t/mean_t);
            circuit.append_op("DEPOLARIZE1", {i}, e);
        }
    }
    return circuit;
}

}   // qontra
