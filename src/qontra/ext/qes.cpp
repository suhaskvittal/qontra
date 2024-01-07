/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include "qontra/ext/qes.h"

#include <vtils/two_level_map.h>

#include <set>

namespace qontra {

std::vector<uint64_t>
get_qubits(const qes::Instruction<>& inst) {
    const std::string name = inst.get_name();

    std::vector<uint64_t> qubits;
    if (name == "h" || name == "cx" || name == "cz" || name == "x" || name == "z"
            || name == "s" || name == "measure" || name == "reset" || name == "liswap")
    {
        for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
            qubits.push_back(inst.get<uint64_t>(i));
        }
    }
    return qubits;
}

std::vector<double>
get_errors(const qes::Instruction<>& inst, const ErrorTable& errors) {
    std::vector<double> error_array;
    auto qubits = get_qubits(inst);
    if (is_2q(inst)) {
        for (size_t i = 0; i < qubits.size(); i += 2) {
            uint64_t q1 = qubits[i],
                     q2 = qubits[i+1];
            double e = vtils::tlm_at(errors.op2q, inst.get_name(), std::make_pair(q1, q2));
            error_array.push_back(e);
        }
    } else {
        for (size_t i = 0; i < qubits.size(); i++) {
            uint64_t q = qubits[i];
            double e;
            if (inst.get_name() == "measure") {
                e = 0.5*(errors.m1w0.at(q) + errors.m0w1.at(q));
            } else {
                e = vtils::tlm_at(errors.op1q, inst.get_name(), q);
            }
            error_array.push_back(e);
        }
    }
    return error_array;
}

std::vector<double>
get_latency(const qes::Instruction<>& inst, const TimeTable& timing) {
    std::vector<double> latency;
    auto qubits = get_qubits(inst);
    if (is_2q(inst)) {
        for (size_t i = 0; i < qubits.size(); i += 2) {
            uint64_t q1 = qubits[i],
                     q2 = qubits[i+1];
            double t = vtils::tlm_at(timing.op2q, inst.get_name(), std::make_pair(q1, q2));
            latency.push_back(t);
        }
    } else {
        for (size_t i = 0; i < qubits.size(); i++) {
            uint64_t q = qubits[i];
            double t = vtils::tlm_at(timing.op1q, inst.get_name(), q);
            latency.push_back(t);
        }
    }
    return latency;
}

}   // qontra
