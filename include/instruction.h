/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef INSTRUCTION_h
#define INSTRUCTION_h

#include <map>
#include <string>

#include <stdint.h>

namespace qontra {

class Instruction {
public:
    enum class type {
        gate,
        measurement,
        event_like
    };

    template <class CONTAINER> Instruction(std::string, CONTAINER args);
    template <class ITERATOR> Instruction(std::string, ITERATOR begin, ITERATOR end);

    Instruction(const Instruction&);
    Instruction(Instruction&&);

    Instruction& operator=(const Instruction&);
};

Instruction::type get_instruction_type(std::string);

typedef std::vector<Instruction> Schedule;

}   // qontra

#include "instruction.inl"

#endif  // INSTRUCTION_H
