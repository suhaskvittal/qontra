/*
 *  author: Suhas Vittal
 *  date:   24 June 2023
 * */

#include "instruction.h"

namespace qontra {

std::vector<uint>
Instruction::get_qubit_operands() {
    // Note: this is for qubit operands that are physically interacted with!
    //
    // Do not include instructions like lckqifmspc which do not interact
    // with qubits directly.

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
        } else if (inst.name == "lckqifmspc") {
            new_inst.operands = std::vector<uint>(inst.operands);
            new_inst.operands[0] = operand_map[inst.operands[0]];
        } else {
            new_inst.operands = inst.operands;
        }
        new_sch.push_back(new_inst);
    }

    return new_sch;
}

}   // qontra
