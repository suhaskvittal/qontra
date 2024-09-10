/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include "qontra/ext/qes.h"

#include <vtils/two_level_map.h>

#include <unordered_set>

namespace qontra {

std::vector<uint64_t>
get_qubits(const qes::Instruction<>& inst) {
    std::vector<uint64_t> qubits;
    if (is_gate(inst)) {
        for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
            qubits.push_back(static_cast<uint64_t>(inst.get<int64_t>(i)));
        }
    }
    return qubits;
}

std::vector<double>
get_errors(const qes::Instruction<>& inst, const ErrorTable& errors) {
    auto qubits = get_qubits(inst);
    if (inst.has_annotation("no_error")) {
        if (is_2q_gate(inst)) {
            return std::vector<double>(qubits.size()/2, 0);
        } else {
            return std::vector<double>(qubits.size(), 0);
        }
    }

    std::vector<double> error_array;
    if (is_2q_gate(inst)) {
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
    auto qubits = get_qubits(inst);
    if (inst.has_annotation("no_tick")) {
        if (is_2q_gate(inst)) {
            return std::vector<double>(qubits.size()/2, 0);
        } else {
            return std::vector<double>(qubits.size(), 0);
        }
    }

    std::vector<double> latency;
    const std::string name = inst.get_name();
    if (is_2q_gate(inst)) {
        if (!timing.op2q.count(name)) {
            std::cerr << "[ get_latency ] no latencies listed for " << name << std::endl;
        }
        for (size_t i = 0; i < qubits.size(); i += 2) {
            uint64_t q1 = qubits[i],
                     q2 = qubits[i+1];
            auto q1_q2 = std::make_pair(q1, q2);
            if (!timing.op2q.at(name).count(q1_q2)) {
                std::cerr << "[ get_latency ] no latency found for "
                    << name << "(" << q1 << ", " << q2 << ")" << std::endl;
                exit(1);
            }
            double t = vtils::tlm_at(timing.op2q, name, q1_q2);
            latency.push_back(t);
        }
    } else {
        if (!timing.op1q.count(name)) {
            std::cerr << "[ get_latency ] no latencies listed for " << name << std::endl;
        }
        for (size_t i = 0; i < qubits.size(); i++) {
            uint64_t q = qubits[i];
            if (!timing.op1q.at(name).count(q)) {
                std::cerr << "[ get_latency ] no latency found for "
                    << name << "(" << q << ")" << std::endl;
                exit(1);
            }
            double t = vtils::tlm_at(timing.op1q, name, q);
            latency.push_back(t);
        }
    }
    return latency;
}

}   // qontra
