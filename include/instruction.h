/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include "defs.h"

#include <algorithm>
#include <string>
#include <vector>

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
};

struct Instruction {
    std::string name;
    std::vector<uint> operands;

    std::vector<uint64_t> exclude_trials;
};

typedef std::vector<Instruction>    schedule_t;

}   // qontra

#endif  // INSTRUCTION_H
