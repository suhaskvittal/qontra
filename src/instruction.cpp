/*
 *  author: Suhas Vittal
 *  date:   24 June 2023
 * */

#include "instruction.h"

struct __asm_inst_t     ASMParserSchedule[4096];
uint32_t                ASMParserScheduleLen = 0;

namespace qontra {

std::string
schedule_to_text(const schedule_t& sch) {
    std::string out;
    for (const auto& inst : sch) {
        out += inst.str() + "\n";
    }
    return out;
}

uint
from_stim_circuit(const stim::Circuit& circ, schedule_t& sch) {
    uint max_qubit = 0;

    uint64_t ectr = 0;  // Event and observable counters.
    uint64_t octr = 0;
    auto cb = [&] (const stim::Operation& op)
    {
        auto gate = op.gate;
        std::string opname(gate->name);
        std::vector<uint> operands;
        // Parse operands:
        //  We need to handle detector and observable instructions
        //  differently, as they are a record lookback (like us!).
        if (opname == "DETECTOR")           operands.push_back(ectr++);
        if (opname == "OBSERVABLE_INCLUDE") operands.push_back(octr++);
        for (auto target : op.target_data.targets) {
            uint32_t x = target.data;
            uint32_t v = x & stim::TARGET_VALUE_MASK;
            operands.push_back(v);
        }
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
            name = "reset";
        } else                              name = "nop";
        if (name == "nop")  return;
        
        sch.push_back({name, operands, {}});
        // Update max_qubit.
        const std::set<std::string> gates{"x", "z", "cx", "h", "reset", "s",
                                        "mrc"};
        if (gates.count(name)) {
            for (uint x : operands) {
                if (x > max_qubit)  max_qubit = x;
            }
        }
    };
    circ.for_each_operation(cb);
    return max_qubit+1; // The qubits were 0-indexed.
}

schedule_t
from_file(std::string fname) {
    FILE* fin = fopen(fname.c_str(), "r");
    asm_yystart(fin);
    asm_yyparse();
    fclose(fin);
    // Convert C-like schedule to C++.
    schedule_t sch;

    for (uint i = 0; i < ASMParserScheduleLen; i++) {
        __asm_inst_t x = ASMParserSchedule[i];
        std::string name(x.name);
        std::vector<uint> operands;
        operands.assign(x.operands.data, x.operands.data + x.operands.size);
        sch.push_back({name, operands, {}});
    }

    return sch;
}

}   // qontra
