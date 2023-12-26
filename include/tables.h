/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef TABLES_h
#define TABLES_h

#include "defs.h"
#include "defs/two_level_map.h"

#include <array>
#include <map>
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

    std::map<uint, fp_t> idling;

    std::map<uint, fp_t> m1w0;
    std::map<uint, fp_t> m0w1;

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
set_all_1q(uint n_qubits, fp_t value, std::map<uint, fp_t>& row) {
    for (uint i = 0; i < n_qubits; i++) {
        row[i] = value;
    }
}

inline void
set_all_2q(uint n_qubits, fp_t value, std::map<std::pair<uint, uint>, fp_t>& row) {
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
    fp_t e_m1w0 = 1e-3; // Readout error (read 0 as 1)
    fp_t e_m0w1 = 1e-3; // Readout error (read 1 as 0)
    fp_t e_idle = 1e-3;
    fp_t t_g1q = 30;    // in nanoseconds.
    fp_t t_g2q = 40;
    fp_t t_ro = 600;
    fp_t t1 = 1000e3;
    fp_t t2 = 500e3;

    ErrorAndTiming& operator*=(fp_t x);
};

inline ErrorAndTiming operator*(ErrorAndTiming et, fp_t x) { et *= x; return et; }
inline ErrorAndTiming operator*(fp_t x, ErrorAndTiming et) { return et * x; }

void populate(uint n_qubits, ErrorTable& errors, TimeTable& timing, const ErrorAndTiming& params);

}   // tables
}   // qontra

#endif // TABLES_h
