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
    "sdg",
    "t",
    "tdg",
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
    "xorfr",    // xorfr   obsno, frno
    // Simulation Instructions
    //
    // Such instructions depend on the simulation being performed.
    // For instance, error enumeration instructions are only used by
    // the error enumerator, but will be treated as NOPs in the a control
    // simulation.
    // 
    // For the error enumerator:
    "cmpx",     // cmpx q1, q2, q3, ... (this is usually an X observable)
    "cmpz",     // cmpz q1, q2, q3, ...
};

const std::set<std::string> IS_QUANTUM {
    "h", "x", "z", "cx", "rz", "s", "sdg", "t", "tdg", "measure", "reset"
};

const std::set<std::string> IS_QUANTUM_LIKE {
    "cmpx", "cmpz"
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

// Annotations and Properties can be defined
// in the program. While there are a set of
// recognized annotations, we leave the use
// of annotations/properties to the user's
// discretion to allow for customization.
//
// Default annotations that are supported:
//      no_tick: ignore instruction latency
//      no_error: ignore instruction error
//      inject_timing_error: inject decoherence and dephasing errors
//  
//  Default properties that are supported:
//      color: 0 (red), 1 (green), 2 (blue) -- for color codes

#define ANNOT_NO_TICK               "no_tick"
#define ANNOT_NO_ERROR              "no_error"
#define ANNOT_INJECT_TIMING_ERROR   "inject_timing_error"

#define PROPERTY_COLOR      "color"
#define PROPERTY_COLOR_RED      0
#define PROPERTY_COLOR_GREEN    1
#define PROPERTY_COLOR_BLUE     2

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

    struct property_value_t {
        int64_t         ival;
        fp_t            fval;
    };

    std::set<std::string> annotations;
    std::map<std::string, property_value_t> properties;

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
