/*
 *  author: Suhas Vittal
 *  date:   5 October 2023
 * */

#ifndef SIM_ENUMERATOR_h
#define SIM_ENUMERATOR_h

#include "defs.h"
#include "experiments.h"
#include "instruction.h"

#include <mpi.h>
#include <stim.h>

#include <fstream>
#include <iostream>
#include <vector>

// The enumerate_errors function (seen below) exhaustively examines the effects of Pauli errors
// on stabilizer circuits. The function only cares about H and CX gates, which can propagate
// errors or switch X errors to Z errors and vice versa. Measurements, resets, and other gates
// are ignored.

namespace qontra {

namespace enumerator {

// Runtime options:
extern bool G_RECORD_WEIGHT_ONE_ERRORS; // Default = false.

enum class error { x, y, z };

// error_record_t = (instruction-number, Instruction, qubit, pauli_error, x_table, z_table)
typedef std::tuple<uint64_t, Instruction, uint64_t, error, stim::simd_bits, stim::simd_bits>
        error_record_t;

void
execute_gate(const Instruction&, stim::simd_bit_table& x_table, stim::simd_bit_table& z_table);

void
inject_error(uint qubit, stim::simd_bit_table& x_table, stim::simd_bit_table& z_table);

}   // enumerator

std::vector<enumerator::error_record_t>
enumerate_errors(const schedule_t& prog);

void
write_recorded_errors_to(std::ostream&, const std::vector<enumerator::error_record_t>&);

}   // qontra

#endif  // SIM_ENUMERATOR_h
