/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"

<<<<<<< HEAD
#include <algorithm>
=======
#include <stim.h>

#include <deque>
>>>>>>> b7fe2b4c19394c47367f32573cb9420d0bd2c7dc
#include <string>

namespace qontra {

const std::vector<std::string> ISA{
    "H",
    "X",
    "Z",
    "CX",
    "S",
    "Mnrc", // Does not record measurement
    "Mrc",  // Records measurement and places it into a buffer
    "R",
    "NOP",          // Essentially a delay operation
    // Control processor instructions.
    "DECODE",       // Tells decoder to decode syndrome.
    "BRDECBUSY",    // Jumps to instruction if decoder is busy (still decoding)
    "FENCEDEC",     // Waits until the decoder finishes.
    "RECORDXOR",    // XORs two entries in the record. Note that operands are
                    // offsets from the end of the array (i.e. 3 means the third
                    // most recent record). First entry is XORd into the second.
    "OBS",          // Computes observable by XORing corresponding results in the
                    // record. Measurement outcome is placed in a buffer.
    "SAVEMEAS",     // The bitstring in the observable buffer is record in a 
                    // probability histogram.
    // Virtual instruction
    "DONE",         // Tells the simulator we are done (virtual instruction)
};

struct Instruction {
    std::string name;
    std::vector<int> operands;

    std::set<uint64_t> exclude_trials;

    // Extra metadata
    bool    is_measuring_x_check;

    std::string str(void) {
        std::string out = name;
        for (uint op : operands)  out += " " + std::to_string(op);
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

#endif  // INSTRUCTION_H
