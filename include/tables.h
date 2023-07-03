/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef TABLES_h
#define TABLES_h

#include "defs.h"

#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace qontra {

class ErrorTable {
public:
    ErrorTable()
        :op1q(), op2q(),
        op2q_leakage_injection(),
        op2q_leakage_transport(),
        op2q_crosstalk(),
        op2q_correlated()
    {}

    //      Entries:    Opname      Affected pair          P(IIII),, ..., P(ZZZZ)
    typedef std::tuple<std::string, std::pair<uint, uint>, std::array<fp_t, 255>> 
        corr_t;

    TwoLevelMap<std::string, uint, fp_t> op1q;
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q;

    // Other errors: not set by default.
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q_leakage_injection;
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q_leakage_transport;
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q_crosstalk;
    TwoLevelMap<std::string, std::pair<uint, uint>, std::vector<corr_t>> op2q_correlated;
};

class TimeTable {  // Units are in nanoseconds.
public:
    TimeTable()
        :op1q(), op2q(), t1(), t2() {}

    TwoLevelMap<std::string, uint, fp_t> op1q;
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q;
    std::map<uint, fp_t> t1;
    std::map<uint, fp_t> t2;
};

namespace tables {

inline void
set_all_1q(uint n_qubits, 
        fp_t value, 
        std::map<uint, fp_t>& row)
{
    for (uint i = 0; i < n_qubits; i++) {
        row[i] = value;
    }
}

inline void
set_all_2q(uint n_qubits,
        fp_t value,
        std::map<std::pair<uint, uint>, fp_t>& row)
{
    for (uint i = 0; i < n_qubits; i++) {
        for (uint j = 0; j < n_qubits; j++) {
            if (i == j) continue;
            auto i_j = std::make_pair(i, j);
            row[i_j] = value;
        }
    }
}

struct ErrorAndTiming {
    fp_t e_g1q = 1e-3;
    fp_t e_g2q = 1e-3;
    fp_t e_ro = 1e-3;
    fp_t t_g1q = 30;    // in nanoseconds.
    fp_t t_g2q = 40;
    fp_t t_ro = 600;
    fp_t t1 = 1000e3;
    fp_t t2 = 500e3;
};

inline void
populate(
        uint n_qubits, 
        ErrorTable& errors,
        TimeTable& timing,
        const ErrorAndTiming& params)
{
    const std::vector<std::string> g1q{"h", "x", "s", "z", "reset"};
    const std::vector<std::string> g2q{"cx"};

    set_all_1q(n_qubits, params.t1, timing.t1);
    set_all_1q(n_qubits, params.t2, timing.t2);
    for (auto g : g1q) {
        set_all_1q(n_qubits, params.e_g1q, errors.op1q[g]);
        set_all_1q(n_qubits, params.t_g1q, timing.op1q[g]);
    }
    // Set measurement characteristics independently.
    set_all_1q(n_qubits, params.e_ro, errors.op1q["m"]);
    set_all_1q(n_qubits, params.t_ro, timing.op1q["m"]);
    for (auto g : g2q) {
        set_all_2q(n_qubits, params.e_g2q, errors.op2q[g]);
        set_all_2q(n_qubits, params.t_g2q, timing.op2q[g]);
    }
}

}   // tables

}   // qontra

#endif // TABLES_h
