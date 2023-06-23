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
    "Mnrc", // does not record measurement
    "Mrc",  // records measurement and places it into a buffer
    "R"
    "DECODE",       // Tells decoder to decode syndrome.
    "BRDECBUSY",    // Jumps to instruction if decoder is busy (still decoding)
    "FENCEDEC",     // Waits until the decoder finishes.
    "NOP",          // essentially a delay operation
    // Virtual instruction
    "DONE",         // Tells the simulator we are done (virtual instruction)
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;

    std::vector<uint64_t> exclude_trials;

    // Extra metadata
    bool    is_measuring_x_check;

    std::string str(void) {
        std::string out = name;
        for (uint op : operands)  out += " " + std::to_string(op);
        return out;
    }
};

typedef std::vector<Instruction>    schedule_t;

#endif  // INSTRUCTION_H
