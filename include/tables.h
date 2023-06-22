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
private:
    fp_t def_1q = 0.0;
    fp_t def_2q = 0.0;
};

class TimeTable {  // Units are in nanoseconds.
public:
    TimeTable()
        :op1q(), op2q(), t1(), t2() {}

    TwoLevelMap<std::string, uint, fp_t> op1q;
    TwoLevelMap<std::string, std::pair<uint, uint>, fp_t> op2q;
    std::map<uint, fp_t> t1;
    std::map<uint, fp_t> t2;
private:
    fp_t def_1q = 30;
    fp_t def_2q = 40;
    fp_t def_t1 = 15000;
    fp_t def_t2 = 7500;

    const fp_t def_ro = 620;
};

} // qontra

#endif // TABLES_h
