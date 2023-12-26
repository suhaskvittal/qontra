/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "tables.h"

#include <limits>

namespace qontra {
namespace tables {

ErrorAndTiming&
ErrorAndTiming::operator*=(fp_t x) {
    this->e_g1q *= x;
    this->e_g2q *= x;
    this->e_m1w0 *= x;
    this->e_m0w1 *= x;
    this->e_idle *= x;
    this->t1 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    this->t2 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    return *this;
}

void populate(uint n_qubits, ErrorTable& errors, TimeTable& timing, const ErrorAndTiming& params) {
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
