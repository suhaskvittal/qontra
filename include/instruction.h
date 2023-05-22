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
    "Z",    // virtual (assuming virtual RZ)
    "CX",
    "S"     // virtual (assuming virtual RZ)
    "M",
    "R"
};

const std::vector<std::string> CP_ISA{
    "INIT",     // Initializes logical qubit
    "EXTRACT",  // Executes syndrome extraction
    "LSMERGE",  // Lattice Surgery Merge
    "LSSPLIT",  // Lattice Surgery Split
    // Physical operations
    "H",
    "X",
    "Z",
    "S",
    "Mnrc",     // Does not record measurement
    "Mrc",      // Records measurement and places it in syndrome buffer.
    "R",
    // Logical operations
    "LH",
    "LX",
    "LZ",
    "LS",
    "LMX",      // qubit, register
    "LMZ",      // qubit, register
    "LR",
    // Compute operations (all integer)
    "ADD",      // register, register, register -- 1 cycle
    "SUB",      // register, register, register -- 1 cycle
    "MUL",      // register, register, register -- 2 cycle pipelined
    "DIV",      // register, register, register -- 4 cycle pipelined
    "NOT",      // register, register           -- 1 cycle
    "AND",      // register, register, register -- 1 cycle
    "OR",       // register, register, register -- 1 cycle
    "XOR",      // register, register, register -- 1 cycle
    // Memory operations
    "LI",       // register, immediate
    "ST",       // register, location
    "LD",       // register, location
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

struct Instruction {   // Physical instruction.
                        // This would be executed by physical operations
                        // on the FTQC.
    std::string name;
    std::vector<uint> operands;
    uint n_ops; // Number of operands. This is not operand.size(): recall that
                // operands is a SIMD-style definition, so we just asking if
                // it is a 1-qubit, 2-qubit, etc. operation.

    std::vector<uint64_t> exclude_trials;
};

}   // qc

namespace cp {

struct Instruction {   // Logical instruction.
                        // This would be executed by Pauli frame updates
                        // or lattice surgery operations.
    std::string name;
    std::vector<uint> operands;
    uint n_ops;

    std::vector<uint64_t> exclude_trials; 
};

}   // cp

#endif  // INSTRUCTION_H
