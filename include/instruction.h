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
    "t",
    "rz",
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
    "any",      // any     obsno, m1, m2, ...
    "xorfr"     // xorfr   obsno, frno
};

const std::set<std::string> IS_QUANTUM_INSTRUCTION {
    "h", "x", "z", "cx", "rz", "t", "s", "measure", "reset"
};

// These operations only have qubits in its operand list.
const std::set<std::string> ONLY_HAS_QUBIT_OPERANDS {
    "h", "x", "z", "s", "t", "cx", "measure", "reset"
};

// These operations have both angles and qubits in the
// instruction list.
const std::set<std::string> INSTRUCTION_USES_ANGLES {
    "rz"
};

const std::set<std::string> IS_2Q_OPERATOR {
    "cx"
};

const std::set<std::string> IS_BR {
    "brifone",
    "brifzero"
};

// These instructions have operands of the form: if-label, event1
const std::set<std::string> IS_BR_TYPE1 {
    "brifone",
    "brifzero"
};

// These instructions are of the form: eventno, measno, measno, ...
const std::set<std::string> IS_DECODE_TYPE1 {
    "event"
};

// These instructions are of the form: obsno, measno, measno, ...
const std::set<std::string> IS_DECODE_TYPE2 {
     "any", "obs"
};

// These instructions are of the form: obsno, frno
const std::set<std::string> IS_DECODE_TYPE3 {
    "xorfr"
};

// These instructions are of the form: frno
const std::set<std::string> IS_DECODE_TYPE4 {
    "decode"
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

    struct operand_block_t {
        std::vector<uint> qubits;
        std::vector<uint> measurements;
        std::vector<uint> events;
        std::vector<uint> observables;
        std::vector<uint> frames;

        std::vector<uint32_t> labels;

        std::vector<fp_t> angles;

        bool operator==(const operand_block_t& other) const {
            return qubits == other.qubits
                    && measurements == other.measurements
                    && events == other.events
                    && observables == other.observables
                    && frames == other.frames
                    && labels == other.labels
                    && angles == other.angles;
        }

        bool operator<(const operand_block_t& other) const {
            if (qubits != other.qubits) return qubits < other.qubits;
            else if (measurements != other.measurements) return measurements < other.measurements;
            else if (events != other.events) return events < other.events;
            else if (observables != other.observables) return observables < other.observables;
            else if (frames != other.frames) return frames < other.frames;
            else if (labels != other.labels) return labels < other.labels;
            else if (angles != other.angles) return angles < other.angles;

            return true;
        }
    } operands;

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
        if (operands.angles.size()) {
            out += "(" + operands_to_str(operands.angles) + ")";
        }
        if (operands.qubits.size()) {
            out += " ";
            out += operands_to_str(operands.qubits);
        }
        if (operands.measurements.size()) {
            out += " ";
            out += operands_to_str(operands.measurements, "M");
        }
        if (operands.events.size()) {
            out += " ";
            out += operands_to_str(operands.events, "e");
        }
        if (operands.observables.size()) {
            out += " ";
            out += operands_to_str(operands.observables, "o");
        }
        if (operands.frames.size()) {
            out += " ";
            out += operands_to_str(operands.frames, "f");
        }
        if (operands.labels.size()) {
            out += " ";
            out += operands_to_str(operands.labels, "addr");
        }
        out += "\n";
        
        return out;
    }

    bool operator<(const Instruction& other) const {
        return name < other.name || (name == other.name && operands < other.operands);
    }

    bool operator==(const Instruction& other) const {
        return name == other.name && operands == other.operands;
    }
private:
    template <typename T>
    std::string operands_to_str(std::vector<T> arr, std::string prefix="") const {
        std::string out;
        bool first = true;
        for (T x : arr) {
            if (!first) out += ", ";
            out += prefix + std::to_string(x);
            first = false;
        }
        return out;
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
stim::Circuit   schedule_to_stim(const schedule_t&, ErrorTable&, TimeTable&, 
                                    fp_t fix_timing_error_as_depolarizing_error=-1);

// The below functions are for parsing ASM (either in a file or in a string).
schedule_t      schedule_from_file(std::string fname);
schedule_t      schedule_from_text(std::string);
schedule_t      schedule_after_parse(void);

}   // qontra

#endif  // INSTRUCTION_H
