/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"

#include <algorithm>
#include <string>
#include <vector>

namespace qontra {

const std::vector<std::string> FTQC_ISA{
    "H",
    "X",
    "Z",
    "CX",
    "S",
    "Mnrc", // does not record measurement
    "Mrc",  // records measurement and places it into the jj buffer
    "R"
};

// The control processor implements a RISCV ISA with quantum extensions.
// Instructions are 64-bit (RISCV will not use upper 32 bits).
//
const std::vector<std::vector<std::string>> CP_ISA{
    // RISV-V Instructions
    // register, register, immediate:
    {"ADDI", "SLTI", "ANDI", "ORI", "XORI", "SLLI", "SRLI", "SRAI", 
        "MULI", "DIVI"},
    // register, immediate
    {"LUI", "AUIPC"},
    // register, register, register
    {"ADD", "SLT", "SLTU", "AND", "OR", "XOR", "SLL", "SRL", "SUB", "SRA",
        "MUL", "DIV"},
    // unconditional jump -- dst register, address
    {"JAL"},
    // jump to base register + offset -- dst register, base register, offset
    {"JALR"}, 
    // branches -- register, register, offset
    {"BEQ", "BNE", "BLT", "BGE"},
    // loads and stores -- src/dst register, base register, address
    {"LDW", "SDW", "LW", "SW", "LH", "SH", "LB", "SB"},
    // Ordering instruction
    {"FENCE"},
    // Exit
    {"EXIT"},
    // Quantum Instructions
    // bookkeeping -- LQ, [PQ] : very large instructions, but are 
    //                              not in the critical path
    //                              necessary for initialization
    {"INITX", "INITZ", "ZOBS", "XOBS"},
    // logical single-qubit instructions: LQ
    {"EXTRACT", "LH", "LX", "LZ", "LS", "LR"},
    // logical multi-qubit instructions: LQ, LQ
    {"LSMERGE", "LSSPLIT"},
    // logical measurements: LQ, register
    {"LMX", "LMZ"},
    // physical single-qubit instructions: PQ
    {"PH", "PX", "PZ", "PS", "PR"},
    // physical multi-qubit instructions: PQ, PQ
    {"PCX"},
    // physical measurements: PQ, register
    {"PM"}
};

namespace qc {
// We define a namespace because control processors
// also have instruction sets

struct Instruction {    // Physical instruction.
                        // This would be executed by physical operations
                        // on the FTQC.
    std::string name;
    std::vector<uint> operands;

    std::vector<uint64_t> exclude_trials;
};

}   // qc

namespace cp {

struct Instruction {    // Logical instruction.
                        // This would be executed by Pauli frame updates
                        // or lattice surgery operations.
    std::string name;
    std::vector<uint> operands;

    std::vector<uint64_t> exclude_trials; 
};

typedef struct {
    cp::Instruction inst;
    bool            with_prev;  // Set if this instruction is in the same
                                // word as the previous instruction in a
                                // schedule.
} vliw_t;

}   // cp

template <class I_t>
using schedule_t = std::vector<I_t>;

}   // qontra

#endif  // INSTRUCTION_H
