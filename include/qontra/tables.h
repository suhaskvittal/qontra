/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef TABLES_h
#define TABLES_h

#include "qontra/defs.h"

#include <vtils/two_level_map.h>

#include <unordered_map>
#include <string>
#include <utility>

namespace qontra {

class ErrorTable {
public:
    ErrorTable() = default;
    ErrorTable(const ErrorTable&) = default;
    ErrorTable(ErrorTable&&) = default;

    ErrorTable& operator=(const ErrorTable&);
    ErrorTable& operator*=(fp_t);

    vtils::TwoLevelMap<std::string, uint64_t, fp_t> op1q;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> op2q;

    std::unordered_map<uint64_t, fp_t> idling;

    std::unordered_map<uint64_t, fp_t> m1w0;
    std::unordered_map<uint64_t, fp_t> m0w1;

    // Other errors: not set by default.
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> 
        op2q_leakage_injection;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> 
        op2q_leakage_transport;
};

class TimeTable {  // Units are in nanoseconds.
public:
    TimeTable() = default;
    TimeTable(const TimeTable&) = default;
    TimeTable(TimeTable&&) = default;

    TimeTable& operator=(const TimeTable&);
    TimeTable& operator*=(fp_t);

    vtils::TwoLevelMap<std::string, uint64_t, fp_t> op1q;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> op2q;
    std::unordered_map<uint64_t, fp_t> t1;
    std::unordered_map<uint64_t, fp_t> t2;
};

ErrorTable operator*(ErrorTable, fp_t);
ErrorTable operator*(fp_t, ErrorTable);

TimeTable operator*(TimeTable, fp_t);
TimeTable operator*(fp_t, TimeTable);

namespace tables {

struct ErrorAndTiming {
    fp_t e_g1q = 1e-3;
    fp_t e_g2q = 1e-3;
    fp_t e_m1w0 = 1e-3; // Readout error (read 0 as 1)
    fp_t e_m0w1 = 1e-3; // Readout error (read 1 as 0)
    fp_t e_idle = 1e-3;
    fp_t t_g1q = 30;    // in nanoseconds.
    fp_t t_g2q = 40;
    fp_t t_ro = 800;
    fp_t t1 = 1000e3;
    fp_t t2 = 500e3;

    ErrorAndTiming& operator*=(fp_t);
};

ErrorAndTiming  operator*(ErrorAndTiming et, fp_t x);
ErrorAndTiming  operator*(fp_t x, ErrorAndTiming et);

void populate(uint64_t n_qubits, ErrorTable&, TimeTable&, const ErrorAndTiming&);

}   // tables
}   // qontra

#include "tables.inl"

#endif // TABLES_h
