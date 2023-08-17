/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"
#include "tables.h"

#include <stim.h>

#include <deque>
#include <string>

#include <stdio.h>

namespace qontra {

const std::vector<std::string> ISA{
    // Quantum Processor instructions
    "h",
    "x",
    "z",
    "cx",
    "s",
    "measure"
    "reset",
    "nop",
    // Control Instructions
    "brifone",    // Branches if measurement is 1 (syntax: brifone LABEL, m-number)
    "brifzero",
    // Decoding Instructions
    "event",    // Creates detection event from args
};

// These instructions have labels in their first operand.
const std::set<std::string> IS_BR_TYPE1 {
    "brifone",
    "brifzero"
};

enum class Annotation {
    // Error annotations
    no_errors,
    inject_timing_errors,
    // Bookkeeping annotations
    round_start,
    flag
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;
    std::vector<Annotation> annotations;

    // Extra metadata (not necessary for general use).
    //
    // Usually, applications will fill this out for their own use.
    // The user should not need to modify this.
    struct {
        uint64_t    owning_check_id = 0;    // The id of the tanner graph vertex
                                            // for which this Instruction is executing
                                            // the check.
        bool        is_for_flag = false;

        std::vector<std::pair<uint64_t, bool>>  operators;  
                                                    // i.e. X1X2X3X4 would be 
                                                    // (1, true), (2, true), etc.
    } metadata;

    std::vector<uint>   get_qubit_operands(void) const;

    std::string str(void) const {
        std::string out = name;
        bool first = true;
        for (uint op : operands) {
            if (first)  out += "\t";
            else        out += ", ";
            out += std::to_string(op);
            first = false;
        }
        return out;
    }

    bool operator<(const Instruction& other) const {
        return name < other.name || (name == other.name && operands < other.operands);
    }

    bool operator==(const Instruction& other) const {
        return name == other.name && operands == other.operands;
    }
};

typedef std::vector<Instruction>    schedule_t;

// Utility functions:
// 
// schedule_to_text converts a schedule to a string (as the name suggests). The
//  string is valid ASM used by the control processor (and thus can be saved to
//  a file).
// schedule_from_stim translates a stim circuit to a schedule. This removes any
//  error signatures (set these in the simulator) and other operations such
//  as coordinate declarations. The output is the max number of qubits in the
//  circuit.

std::string     schedule_to_text(const schedule_t&);

schedule_t      schedule_from_stim(const stim::Circuit&);
schedule_t      relabel_operands(const schedule_t&);
schedule_t      schedule_from_file(std::string fname);
schedule_t      schedule_from_text(std::string);

schedule_t      schedule_after_parse(void);

}   // qontra

#endif  // INSTRUCTION_H
