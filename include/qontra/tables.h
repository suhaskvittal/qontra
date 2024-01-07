/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef TABLES_h
#define TABLES_h

#include "qontra/defs.h"

#include <vtils/two_level_map.h>

#include <map>
#include <string>
#include <utility>

namespace qontra {

class ErrorTable {
public:
    ErrorTable();
    ErrorTable(const ErrorTable&);
    ErrorTable(ErrorTable&&);

    ErrorTable& operator=(const ErrorTable&);

    vtils::TwoLevelMap<std::string, uint64_t, fp_t> op1q;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> op2q;

    std::map<uint64_t, fp_t> idling;

    std::map<uint64_t, fp_t> m1w0;
    std::map<uint64_t, fp_t> m0w1;

    // Other errors: not set by default.
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> 
        op2q_leakage_injection;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> 
        op2q_leakage_transport;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> 
        op2q_crosstalk;
};

class TimeTable {  // Units are in nanoseconds.
public:
    TimeTable();
    TimeTable(const TimeTable&);
    TimeTable(TimeTable&&);

    TimeTable& operator=(const TimeTable&);

    vtils::TwoLevelMap<std::string, uint64_t, fp_t> op1q;
    vtils::TwoLevelMap<std::string, std::pair<uint64_t, uint64_t>, fp_t> op2q;
    std::map<uint64_t, fp_t> t1;
    std::map<uint64_t, fp_t> t2;
};

namespace tables {

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

ErrorAndTiming  operator*(ErrorAndTiming et, fp_t x);
ErrorAndTiming  operator*(fp_t x, ErrorAndTiming et);

void populate(uint64_t n_qubits, ErrorTable&, TimeTable&, const ErrorAndTiming&);

}   // tables
}   // qontra

#include "tables.inl"

#endif // TABLES_h
