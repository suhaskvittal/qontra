/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"

#include <string>
#include <vector>

const std::vector<std::string> FTQC_ISA{
    "H",
    "X",
    "Z",
    "CX",
    "S",
    "Mnrc", // does not record measurement
    "Mrc",  // records measurement and places it into the syndrome buffer
    "R"
};

const std::vector<std::string> CP_ISA{
    "INIT",     // logical-qubit, [qubits]  -- initializes logical qubit,
                //                              even qubits = parity
                //                              odd qubits = data
    "EXTRACT",  // logical-qubit            -- performs syndrome extraction
    "LSMERGE",  // qubit, qubit             -- performs lattice surgery merge
    "LSSPLIT",  // qubit, qubit             -- performs lattice surgery split
    // Physical operations
    "H",
    "X",
    "Z",
    "S",
    "CX",
    "M",
    "R",
    // Logical operations (SIMD)
    "LH",       // qubit
    "LX",       // qubit
    "LZ",       // qubit
    "LS",       // qubit
    "LMX",      // qubit, register
    "LMZ",      // qubit, register
    // Classical operations (RISC-Style)
    // Three types:
    //  R2: r2 = op r1
    //  R3: r3 = r2 op r1
    //  RI: register + immediate
    //  BO: base + offset
    //  
    //  Specifiying a register $0 corresponds to the zero register.
    //  Standard data type is integer.
    "ADD",      // register, register, register
    "SUB",      // register, register, register
    "MUL",      // register, register, register
    "DIV",      // register, register, register
    "NOT",      // register, register          
    "AND",      // register, register, register
    "OR",       // register, register, register
    "XOR",      // register, register, register
    "LI",       // register, immediate
    "ST",       // register[src], location, register[offset]
    "LD",       // register[dst], location, register[offset]
    // Pseudo operations (cannot be executed on a real system, but
    // are useful in simulation)
    "PEEKX",    // Check X Pauli frame with ideal decoder (cannot err), 
                //  place in register
    "PEEKZ"     // Check Z Pauli frame with ideal decoder (cannot err),
                //  place in register
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

    enum class Type {r2, r3, bo, ri};
    Type instruction_type;

    std::vector<uint64_t> exclude_trials; 
};

}   // cp

#endif  // INSTRUCTION_H
