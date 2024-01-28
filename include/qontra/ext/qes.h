/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#ifndef QONTRA_EXT_QES_h
#define QONTRA_EXT_QES_h

#include "qontra/isa.h"
#include "qontra/tables.h"

#include <qes.h>

#include <vector>

namespace qontra {

std::vector<uint64_t>   get_qubits(const qes::Instruction<>&);
std::vector<double>     get_errors(const qes::Instruction<>&, const ErrorTable&);
std::vector<double>     get_latency(const qes::Instruction<>&, const TimeTable&);

double  get_max_latency(const qes::Instruction<>&, const TimeTable&);

isa_data_t  isa_get(const qes::Instruction<>&);

bool    is_gate(const qes::Instruction<>&);
bool    is_2q_gate(const qes::Instruction<>&);

bool    is_instantaneous(const qes::Instruction<>&);

size_t  get_number_of_qubits(const qes::Program<>&);
size_t  get_number_of_qubits(const qes::Program<>&, std::vector<uint64_t>& qubits);

}   // qontra

#include "qes.inl"

#endif  // QONTRA_EXT_QES_h
