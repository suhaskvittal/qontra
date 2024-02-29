/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "qontra/tables.h"

#include <limits>
#include <vector>

namespace qontra {

template <class K> inline void
inplace_mul(std::map<K, fp_t>& m, fp_t x) {
    for (auto& [k, p] : m) {
        p *= x;
    }
}

template <class K> inline void
inplace_mul(vtils::TwoLevelMap<std::string, K, fp_t>& m, fp_t x) {
    for (auto& [s, _m] : m) {
        inplace_mul(_m, x);
    }
}

ErrorTable&
ErrorTable::operator*=(fp_t x) {
    inplace_mul(op1q, x);
    inplace_mul(op2q, x);
    inplace_mul(idling, x);
    inplace_mul(m1w0, x);
    inplace_mul(m0w1, x);
    inplace_mul(op2q_leakage_injection, x);
    inplace_mul(op2q_leakage_transport, x);
    return *this;
}

TimeTable&
TimeTable::operator*=(fp_t x) {
    inplace_mul(t1, x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x);
    inplace_mul(t2, x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x);
    return *this;
}

namespace tables {

ErrorAndTiming&
ErrorAndTiming::operator*=(fp_t x) {
    e_g1q *= x;
    e_g2q *= x;
    e_m1w0 *= x;
    e_m0w1 *= x;
    e_idle *= x;
    t1 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    t2 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    return *this;
}

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
