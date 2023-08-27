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
#include <math.h>

namespace qontra {

const std::vector<std::string> ISA{
    // Quantum Processor instructions
    "h",
    "x",
    "z",
    "cx",
    "s",
    "measure",
    "reset",
    "nop",
    // Control Instructions
    "brifone",  // brifone LABEL, event
    "brifzero",
    // Decoding Instructions
    "decode",   // decode  frno
    "event",    // event   eventno, m1, m2, m3, ...
    "obs",      // obs     obsno, m1, m2, ...
    "xorfr"     // xorfr   obsno, frno
};

const std::set<std::string> IS_QUANTUM_INSTRUCTION {
    "h", "x", "z", "cx", "s", "measure", "reset"
};

const std::set<std::string> ONLY_HAS_QUBIT_OPERANDS {
    "h", "x", "z", "s", "cx", "measure", "reset"
};

const std::set<std::string> IS_2Q_OPERATOR {
    "cx"
};

// These instructions have labels in their first operand.
const std::set<std::string> IS_BR_TYPE1 {
    "brifone",
    "brifzero"
};

// Each instruction can have annotations that indicate
// properties of the program state. 
//
// no_error: do not inject errors for this instruction
// no_tick: ignore decohering/dephasing effects of this instruction
// inject_timing_errors: inject decoherence/dephasing errors at this instruction.
//
// round_start: round starts at this instruction
// flag: this is a flag measurement

enum class Annotation {
    // Error annotations
    no_error,
    no_tick,
    inject_timing_error,
    // Bookkeeping annotations
    round_start,
    flag
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;
    std::set<Annotation> annotations;

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
            if (first)  out += " ";
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

uint            get_number_of_qubits(const schedule_t&);

std::string     schedule_to_text(const schedule_t&);
stim::Circuit   schedule_to_stim(const schedule_t&, ErrorTable&, TimeTable&);

// The below functions are for parsing ASM (either in a file or in a string).
schedule_t      schedule_from_file(std::string fname);
schedule_t      schedule_from_text(std::string);
schedule_t      schedule_after_parse(void);

}   // qontra

#endif  // INSTRUCTION_H
