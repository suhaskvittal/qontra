/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "qontra/tables.h"

#include <limits>
#include <vector>

namespace qontra {

ErrorTable::ErrorTable()
    :op1q(),
    op2q(),
    idling(),
    m1w0(),
    m0w1(),
    op2q_leakage_injection(),
    op2q_leakage_transport(),
    op2q_crosstalk()
{}

ErrorTable::ErrorTable(const ErrorTable& other)
    :op1q(other.op1q),
    op2q(other.op2q),
    idling(other.idling),
    m1w0(other.m1w0),
    m0w1(other.m0w1),
    op2q_leakage_injection(other.op2q_leakage_injection),
    op2q_leakage_transport(other.op2q_leakage_injection),
    op2q_crosstalk(other.op2q_leakage_injection)
{}

ErrorTable::ErrorTable(ErrorTable&& other)
    :op1q(std::move(other.op1q)),
    op2q(std::move(other.op2q)),
    idling(std::move(other.idling)),
    m1w0(std::move(other.m1w0)),
    m0w1(std::move(other.m0w1)),
    op2q_leakage_injection(std::move(other.op2q_leakage_injection)),
    op2q_leakage_transport(std::move(other.op2q_leakage_injection)),
    op2q_crosstalk(std::move(other.op2q_leakage_injection))
{}

TimeTable::TimeTable()
    :op1q(),
    op2q(),
    t1(),
    t2()
{}

TimeTable::TimeTable(const TimeTable& other)
    :op1q(other.op1q),
    op2q(other.op2q),
    t1(other.t1),
    t2(other.t2)
{}

TimeTable::TimeTable(TimeTable&& other)
    :op1q(std::move(other.op1q)),
    op2q(std::move(other.op2q)),
    t1(std::move(other.t1)),
    t2(std::move(other.t2))
{}

namespace tables {

inline void
set_all_1q(uint64_t n_qubits, fp_t value, std::map<uint64_t, fp_t>& row) {
    for (uint64_t i = 0; i < n_qubits; i++) {
        row[i] = value;
    }
}

inline void
set_all_2q(uint64_t n_qubits, fp_t value, std::map<std::pair<uint64_t, uint64_t>, fp_t>& row) {
    for (uint64_t i = 0; i < n_qubits; i++) {
        for (uint64_t j = 0; j < n_qubits; j++) {
            if (i == j) continue;
            auto i_j = std::make_pair(i, j);
            row[i_j] = value;
        }
    }
}

void populate(uint64_t n_qubits, ErrorTable& errors, TimeTable& timing, const ErrorAndTiming& params) {
    // Assume z, rz, s, sdg, t, tdg are implemented via a virtual RZ.
    const std::vector<std::string> g1q{"h", "x", "rx", "ry", "reset"};
    const std::vector<std::string> g2q{"cx", "liswap"};

    set_all_1q(n_qubits, params.t1, timing.t1);
    set_all_1q(n_qubits, params.t2, timing.t2);
    for (std::string g : g1q) {
        set_all_1q(n_qubits, params.e_g1q, errors.op1q[g]);
        set_all_1q(n_qubits, params.t_g1q, timing.op1q[g]);
    }
    set_all_1q(n_qubits, params.e_idle, errors.idling);
    // Set measurement characteristics independently.
    set_all_1q(n_qubits, params.e_m1w0, errors.m1w0);
    set_all_1q(n_qubits, params.e_m0w1, errors.m0w1);
    set_all_1q(n_qubits, params.t_ro, timing.op1q["measure"]);
    for (std::string g : g2q) {
        set_all_2q(n_qubits, params.e_g2q, errors.op2q[g]);
        set_all_2q(n_qubits, params.t_g2q, timing.op2q[g]);
    }
}

}   // tables
}   // qontra
