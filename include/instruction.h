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
    // Qubit locking instructions. If a qubit is locked, any operation that
    // interacts with it physically is treated as a nop.
    "lockq",
    "unlockq",
    // Jumps + Branches + Conditions
    "jmp"           // Jumps to address unconditionally
    "brdb",         // Jumps to address if decoder is busy (still decoding)
    // Event speculation instructions
    "brifmspc",     // Jumps to address if all measurements in list have been
                    // speculated.
    "brnifmspc",    // Same as brifmspc, but if not all have been speculated.
    "mspcflush",    // Clears all measurement speculations for the provided
                    // operands.
    // Ordering instructions
    "dfence",       // Waits until decoder finishes.
    // Tracking instructions
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
    "hshift",       // Shifts event history by the specified operand (oldest "k"
                    // events are removed).
    // Virtual instruction
    "done",         // Tells the simulator we are done
    "ffstart",      // Tells the simulator to start fast-forwarding
    "ffend"         // Tells the simulator to stop fast-forwarding
};

const std::set<std::string> ONLY_HAS_QUBIT_OPERANDS{
    "h",
    "x",
    "z",
    "cx",
    "s",
    "mnrc",
    "reset",
    "nop"
};

// These instructions should wait until all qubit operations
// have finished. This is typically if they rely on measurements
// finishing.
const std::set<std::string> IS_FENCE{
    "event",
    "obs"
};

// These instructions have nothing to be done during the QEX stage.
const std::set<std::string> IS_NOP_LIKE{
    "ffstart",
    "ffend"
};

const std::set<std::string> ARE_JMP_OR_BR{
    "jmp",
    "brdb",
    "brifmspc",
    "brnifmspc",
};

const std::set<std::string> IS_LOCKING{
    "lockq",
    "unlockq",
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;

    // Extra metadata (not necessary for general use).
    //
    // Usually, applications will fill this out for their own use.
    // The user should not need to modify this.
    struct {
        uint64_t    owning_check_id = 0;    // The id of the tanner graph vertex
                                            // for which this Instruction is executing
                                            // the check.
        bool        is_for_flag = false;
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
// fast_convert_to_stim translates a schedule to a stim circuit and injects errors
//  according to the error and timing information provided. Note that this is approximate
//  and assumes that each instruction is a discrete time step. For more accurate
//  error models, use the ControlSimulator's build_error_model function. Note that the
//  ControlSimulator is also much slower.

std::string     schedule_to_text(const schedule_t&);

schedule_t      schedule_from_stim(const stim::Circuit&);
schedule_t      relabel_operands(const schedule_t&);
schedule_t      schedule_from_file(std::string fname);
schedule_t      schedule_from_text(std::string);

schedule_t      schedule_after_parse(void);

stim::Circuit   fast_convert_to_stim(const schedule_t&, ErrorTable&, TimeTable&);

}   // qontra

#endif  // INSTRUCTION_H
