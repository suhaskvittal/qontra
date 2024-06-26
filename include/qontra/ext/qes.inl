/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <algorithm>
#include <set>

namespace qontra {

inline double
get_max_latency(const qes::Instruction<>& inst, const TimeTable& timing) {
    std::vector<double> latencies = get_latency(inst, timing);
    return *std::max_element(latencies.begin(), latencies.end());
}

inline isa_data_t
isa_get(const qes::Instruction<>& inst) {
    return isa_get(inst.get_name());
}

inline bool
is_gate(const qes::Instruction<>& inst) {
    INSTRUCTION_TYPE t = isa_get(inst).inst_type;
    return t == INSTRUCTION_TYPE::QUANTUM_G1Q || t == INSTRUCTION_TYPE::QUANTUM_G2Q;
}

inline bool
is_2q_gate(const qes::Instruction<>& inst) {
    return isa_get(inst).inst_type == INSTRUCTION_TYPE::QUANTUM_G2Q;
}

inline bool
is_instantaneous(const qes::Instruction<>& inst) {
    return !is_gate(inst);
}

inline size_t
get_number_of_qubits(const qes::Program<>& program) {
    std::vector<uint64_t> tmp;
    return get_number_of_qubits(program, tmp);
}

inline size_t
get_number_of_qubits(const qes::Program<>& program, std::vector<uint64_t>& qubits) {
    uint64_t max_qubit = 0;
    std::set<uint64_t> visited;
    for (const auto& inst : program) {
        auto _qubits = get_qubits(inst);
        for (uint64_t q : _qubits) {
            if (!visited.count(q)) {
                qubits.push_back(q);
                visited.insert(q);
            }
        }
        if (_qubits.empty()) continue;
        uint64_t mq = *std::max_element(_qubits.begin(), _qubits.end());
        max_qubit = std::max(max_qubit, mq);
    }
    return max_qubit+1;
}

}   // qontra
