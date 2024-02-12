/*
 *  author: Suhas Vittal
 *  date:   9 February 2024
 * */

#include <qontra/ext/qes.h>

#include <vtils/set_algebra.h>

#include <iostream>
#include <map>
#include <set>

using namespace qontra;

int main(int argc, char* argv[]) {
    std::string qes_file(argv[1]);

    qes::Program<> ext = qes::from_file(qes_file);

    std::set<int64_t> tracked_set;
    std::set<int64_t> obs_qubit_set;
    for (size_t i = 2; i < argc; i++) {
        tracked_set.insert(atoi(argv[i]));
    }
    int64_t old_mctr = 0;
    std::map<int64_t, int64_t> old_mctr_to_qubit;
    // First, identify all qubits that interact with any qubit in the tracked_set.
    // These qubits are added to the tracked_set.
    //
    // We will also track qubits that are measured for the observables, but
    // these are distinct from the tracked_set.
    for (const qes::Instruction<>& inst : ext) {
        std::string name = inst.get_name();
        if (name == "cx") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i += 2) {
                int64_t q1 = inst.get<int64_t>(i),
                        q2 = inst.get<int64_t>(i+1);
                if (tracked_set.count(q2)) tracked_set.insert(q1);
            }
        } else if (name == "measure") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
                int64_t q = inst.get<int64_t>(i);
                old_mctr_to_qubit[old_mctr++] = q;
            }
        } else if (name == "obs") {
            for (size_t i = 1; i < inst.get_number_of_operands(); i++) {
                int64_t m = inst.get<int64_t>(i);
                obs_qubit_set.insert(old_mctr_to_qubit[m]);
            }
        }
    }
    // Now, return a program with operations containing qubits in the tracked
    // set.
    old_mctr = 0;
    int64_t new_mctr = 0;
    std::map<int64_t, int64_t> old_to_new_mctr_map;

    qes::Program<> taint_program;
    for (const qes::Instruction<>& inst : ext) {
        std::string name = inst.get_name();
        std::vector<int64_t> operands;
        if (name == "h" || name == "reset" || name == "measure") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
                int64_t q = inst.get<int64_t>(i);
                if (tracked_set.count(q) || (name == "measure" && obs_qubit_set.count(q))) {
                    operands.push_back(q);
                    if (name == "measure") {
                        old_to_new_mctr_map[old_mctr] = new_mctr++;
                    }
                }
                old_mctr += (name == "measure");
            }
        } else if (name == "cx") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i += 2) {
                int64_t q1 = inst.get<int64_t>(i),
                        q2 = inst.get<int64_t>(i+1);
                if (tracked_set.count(q1) && tracked_set.count(q2)) {
                    if (tracked_set.count(q1) ^ tracked_set.count(q2)) {
                        std::cerr << "[ warning ] CX(" << q1 << ", " << q2 << ") has one qubit untracked: "
                            << q1 << "(" << tracked_set.count(q1) << "), "
                            << q2 << "(" << tracked_set.count(q2) << ")" << std::endl;
                    }
                    operands.push_back(q1);
                    operands.push_back(q2);
                }
            }
        } else {
            // event or obs.
            operands.push_back(inst.get<int64_t>(0));
            for (size_t i = 1; i < inst.get_number_of_operands(); i++) {
                int64_t m = inst.get<int64_t>(i);
                if (old_to_new_mctr_map.count(m)) {
                    operands.push_back(old_to_new_mctr_map[m]);
                } else {
                    operands.clear();
                    goto add_instruction;
                }
            }
        }
add_instruction:
        if (!operands.empty()) {
            taint_program.emplace_back(name, operands);
        }
    }
    std::cout << taint_program << std::endl;
}
