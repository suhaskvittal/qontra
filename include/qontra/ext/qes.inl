/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <set>

namespace qontra {

inline bool
is_gate(const qes::Instruction<>& inst) {
    const std::string name = inst.get_name();
    return name == "h" || name == "x" || name == "z" || name == "cx"
        || name == "cz" || name == "s" || name == "measure" 
        || name == "reset" || name == "liswap";
}

inline bool
is_2q(const qes::Instruction<>& inst) {
    const std::string name = inst.get_name();
    return name == "cx" || name == "cz" || name == "liswap";
}

inline bool
is_instantaneous(const qes::Instruction<>& inst) {
    const std::string name = inst.get_name();
    return name == "event" || name == "obs" || name == "mshift";
}

inline bool
error_goes_before_op(const qes::Instruction<>& inst) {
    return inst.get_name() == "measure";
}

inline bool
apply_x_error_instead_of_depolarizing(const qes::Instruction<>& inst) {
    return inst.get_name() == "measure" || inst.get_name() == "reset";
}

inline size_t
get_number_of_qubits(const qes::Program<>& program) {
    std::vector<uint64_t> tmp;
    return get_number_of_qubits(program, tmp);
}

inline size_t
get_number_of_qubits(const qes::Program<>& program, std::vector<uint64_t>& qubits) {
    std::set<uint64_t> qubits_in_program;
    for (const auto& inst : program) {
        auto _qubits = get_qubits(inst);
        qubits_in_program.insert(_qubits.begin(), _qubits.end());
    }
    qubits = std::vector<uint64_t>(qubits_in_program.begin(), qubits_in_program.end());
    return qubits_in_program.size();
}

}   // qontra
