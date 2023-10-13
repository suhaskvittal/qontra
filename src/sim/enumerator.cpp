/*
 *  author: Suhas Vittal
 *  date:   5 October 2023
 * */

#include "sim/enumerator.h"

namespace qontra {
namespace enumerator {

bool G_RECORD_WEIGHT_ONE_ERRORS = false;

void
execute_gate(const Instruction& inst, stim::simd_bit_table& x_table, stim::simd_bit_table& z_table) {
    std::vector<uint> qubits(inst.operands.qubits);
    // We really only care about H and CX, and the cmpz and cmpx instructions.
    if (inst.name == "h") {
        for (uint q : qubits) {
            x_table[q].swap_with(z_table[q]);
        }
    } else if (inst.name == "cx") {
        for (uint i = 0; i < qubits.size(); i += 2) {
            uint q0 = qubits[i], q1 = qubits[i+1];

            x_table[q1] ^= x_table[q0];
            z_table[q0] ^= z_table[q1];
        }
    } else if (inst.name == "cmpx") {
        for (uint q = 0; q < x_table.num_major_bits_padded(); q++) {
            if (!qubits.empty() && std::find(qubits.begin(), qubits.end(), q) == qubits.end()) {
                x_table[q].clear();
            }
        }
        z_table.clear();
    } else if (inst.name == "cmpz") {
        for (uint q = 0; q < z_table.num_major_bits_padded(); q++) {
            if (!qubits.empty() && std::find(qubits.begin(), qubits.end(), q) == qubits.end()) {
                z_table[q].clear();
            }
        }
        x_table.clear();
    }
}

void
inject_error(uint q, stim::simd_bit_table& x_table, stim::simd_bit_table& z_table) {
    x_table[q][0] = 1;
    x_table[q][1] = 1;
    x_table[q][2] = 0;

    z_table[q][0] = 0;
    z_table[q][1] = 1;
    z_table[q][2] = 1;
}

}   // enumerator

using namespace enumerator;
using namespace experiments;

std::vector<error_record_t>
enumerate_errors(const schedule_t& prog) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }

    const uint n = get_number_of_qubits(prog);
    stim::simd_bit_table x_table(n, 3);
    stim::simd_bit_table z_table(n, 3);

    std::vector<error_record_t> record_list;

    const error error_list[] = { error::x, error::y, error::z };

    uint64_t k = 0;
    for (uint i = 0; i < prog.size(); i++) {
        Instruction inst = prog[i];
        if (inst.name != "h" && inst.name != "cx") {
            continue;
        }
        auto qubits = inst.operands.qubits;
        for (uint ii = 0; ii < qubits.size(); ii++) {
            x_table.clear();
            z_table.clear();
            if ((k++) % world_size != world_rank) {
                continue;
            }
            uint q = qubits[ii];
            inject_error(qubits[ii], x_table, z_table);
            for (uint j = i+1; j < prog.size(); j++) {
                execute_gate(prog[j], x_table, z_table);
            }
            // Now the program should be done -- check the x and z tables for errors.
            auto x_tr = x_table.transposed();
            auto z_tr = z_table.transposed();
            for (uint j = 0; j < 3; j++) {
                stim::simd_bits_range_ref x_ref = x_tr[j];
                stim::simd_bits_range_ref z_ref = z_tr[j];
                
                if (((x_ref.popcnt() + G_RECORD_WEIGHT_ONE_ERRORS) > 1)
                    || (z_ref.popcnt() + G_RECORD_WEIGHT_ONE_ERRORS) > 1)
                {
                    stim::simd_bits x_cpy(x_ref), z_cpy(z_ref);
                    error_record_t rec = std::make_tuple(i, inst, q, error_list[j], x_cpy, z_cpy);
                    record_list.push_back(rec);
                }
            }
        }
    }
    // Collect results.
    std::vector<error_record_t> final_record_list;
    if (G_USE_MPI && world_size > 1) {
        for (int i = 0; i < world_size; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            uint n_recs = record_list.size();
            MPI_Bcast(&n_recs, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
            for (uint j = 0; j < n_recs; j++) {
                error_record_t& rec = record_list[j];

                uint64_t inst_no = std::get<0>(rec);
                MPI_Bcast(&inst_no, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
                Instruction inst = prog[inst_no];

                uint64_t q = std::get<2>(rec);
                MPI_Bcast(&q, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);

                error e = std::get<3>(rec);
                uint ei;    // We need to transmit an integer -- cannot transmit an enum.
                if (e == error::x)      ei = 0;
                else if (e == error::y) ei = 1;
                else                    ei = 2;
                MPI_Bcast(&ei, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
                e = error_list[ei];

                stim::simd_bits x(std::get<4>(rec));
                stim::simd_bits z(std::get<5>(rec));
                MPI_Bcast(x.u64, x.num_u64_padded(), MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
                MPI_Bcast(z.u64, z.num_u64_padded(), MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);

                error_record_t new_rec = std::make_tuple(inst_no, inst, q, e, x, z);
                final_record_list.push_back(new_rec);
            }
        }
    } else {
        final_record_list = record_list;
    }
    return final_record_list;
}

void
write_recorded_errors_to(std::ostream& out, const std::vector<error_record_t>& record_list) {
    for (const error_record_t& rec : record_list) {
        uint64_t inst_no = std::get<0>(rec);
        Instruction inst = std::get<1>(rec);
        uint64_t q = std::get<2>(rec);
        error e = std::get<3>(rec);
        if (e == error::y)  continue;
        stim::simd_bits x = std::get<4>(rec);
        stim::simd_bits z = std::get<5>(rec);

        typedef std::pair<uint, error> qubit_error_t;
        std::vector<qubit_error_t> error_list;

        uint64_t n = x.num_bits_padded();   // Good enough approximation.
        // Cmopute the list of errors that have occurred due to this error.
        for (uint64_t i = 0; i < n; i++) {
            if (x[i] & z[i])    error_list.push_back(std::make_pair(i, error::y));
            else if (x[i])      error_list.push_back(std::make_pair(i, error::x));
            else if (z[i])      error_list.push_back(std::make_pair(i, error::z));
        }
        // For ease, create a map from error to string.
        std::map<error, std::string> error_str_map;
        error_str_map[error::x] = "x";
        error_str_map[error::y] = "y";
        error_str_map[error::z] = "z";
        // Write error and its outcomes to the ostream.
        out << "[ pc = " << std::hex << inst_no << std::dec << " ] instruction: " << inst.name;
        if (inst.name == "cx") {
            for (uint i = 0; i < inst.operands.qubits.size(); i += 2) {
                uint r1 = inst.operands.qubits[i];
                uint r2 = inst.operands.qubits[i+1];
                if (r1 == q || r2 == q) {
                    out << " " << r1 << " " << r2;
                } else {
                    out << " _ _";
                }
            }
        } else {
            for (uint r : inst.operands.qubits) {
                if (r == q) out << " " << q;
                else        out << " _";
            }
        }
        out << "\n\t" << error_str_map[e] << q << "\t";
        for (auto& qe : error_list) {
            uint64_t r = qe.first;
            error f = qe.second;
            out << " " << error_str_map[f] << r;
        }
        out << "\n";
    }
}

}   // qontra
