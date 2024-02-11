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
    std::string protean_folder(argv[1]);

    std::string ext_file = protean_folder + "/ext.qes";
    qes::Program<> ext = qes::from_file(ext_file);
    // Propagate Paulis through circuit.
    std::set<uint64_t> zero_set;
    std::map<uint64_t, std::set<uint64_t>> x_prop_map, z_prop_map;
    
    const uint64_t n = get_number_of_qubits(ext);
    for (uint64_t q = 0; q < n; q++) {
        x_prop_map[q] = {q};
        z_prop_map[q] = {q};
    }

    size_t pc = 0;
    for (const qes::Instruction<>& inst : ext) {
        std::string name = inst.get_name();
        if (name == "reset") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
                uint64_t q = static_cast<uint64_t>(inst.get<int64_t>(i));
                zero_set += q;
                x_prop_map.erase(q);
                z_prop_map.erase(q);
            }
        } else if (name == "h") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
                uint64_t q = static_cast<uint64_t>(inst.get<int64_t>(i));
                zero_set -= q;
            }
        } else if (name == "cx") {
            for (size_t i = 0; i < inst.get_number_of_operands(); i += 2) {
                uint64_t q1 = static_cast<uint64_t>(inst.get<int64_t>(i)),
                         q2 = static_cast<uint64_t>(inst.get<int64_t>(i+1));
                std::cout << "cnot between " << q1 << " and " << q2 << ":\n";
                std::cout << "\tbefore: ";
                for (uint64_t x : x_prop_map[q1]) std::cout << "x" << x;
                std::cout << " ^ ";
                for (uint64_t x : x_prop_map[q2]) std::cout << "x" << x;
                std::cout << "\n\tafter: ";
                for (uint64_t x : x_prop_map[q1]) std::cout << "x" << x;
                std::cout << " , ";
                for (uint64_t x : x_prop_map[q2]) std::cout << "x" << x;
                std::cout << "\n";
                x_prop_map[q1] ^= x_prop_map[q2];
                z_prop_map[q2] ^= z_prop_map[q1];
            }
        } else if (name == "measure") {
            std::cout << "measurement @ pc = " << pc << ":\n";
            for (size_t i = 0; i < inst.get_number_of_operands(); i++) {
                uint64_t q = static_cast<uint64_t>(inst.get<int64_t>(i));
                std::cout << "\t" << q << " --> ";
                for (uint64_t x : x_prop_map[q]) std::cout << "x" << x;
                std::cout << " ";
                for (uint64_t z : z_prop_map[q]) std::cout << "z" << z;
                std::cout << "\n";
            }
        }
        pc++;
    }
}
