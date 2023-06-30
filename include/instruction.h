/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"

#include <stim.h>

#include <deque>
#include <string>

#include <stdio.h>

namespace qontra {

const std::vector<std::string> ISA{
    "h",
    "x",
    "z",
    "cx",
    "s",
    "mnrc", // Does not record measurement
    "mrc",  // Records measurement and places it into a buffer
    "reset",
    "nop",          // Essentially a delay operation
    // Control processor instructions.
    "decode",       // Tells decoder to decode syndrome and places the resulting
                    // Pauli frame at the offset specified by the operand.
    "brdb",         // Jumps to instruction if decoder is busy (still decoding)
    "dfence",       // Waits until decoder finishes.
    "event",        // Creates a detection event and writes to an event history
                    // buffer. The first operand indicates the event index (0
                    // to 4095), and the remaining operands' (which are lookback 
                    // indices into the measurement record) corresponding entries 
                    // are XORd to create the event.
    "obs",          // Computes observable by XORing operands (after the first) 
                    // in the record. Observable is placed in a buffer at the
                    // location specified by the first operand.
    "xorfr",        // Xors the Pauli frame into the specified obs buffer location.
    "savem",        // The bitstring in the observable buffer is record in a 
                    // probability histogram.
    // Virtual instruction
    "done",         // Tells the simulator we are done (virtual instruction)
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;

    std::set<uint64_t> exclude_trials;

    // Extra metadata
    bool    is_measuring_x_check;

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
// schedule_to_text converts a schedule to a string (as the name suggests).
// from_stim_circuit translates a stim circuit to a schedule. This removes any
//  error signatures (set these in the simulator) and other operations such
//  as coordinate declarations. The output is the max number of qubits in the
//  circuit.

std::string     schedule_to_text(const schedule_t&);
uint            from_stim_circuit(const stim::Circuit&, schedule_t&);
schedule_t      from_file(std::string fname);

}   // qontra

#endif  // INSTRUCTION_H
