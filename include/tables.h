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

template <typename T, typename U, typename V>
using map2d<T, U, V> = std::map<T, std::map<U, V>>;

template <typename T, typename U, typename V>
void put(map2d<T, U, V>& x, T i, U j, V k) {
    if (!x.count(i)) {
        x[i] = std::map<U, V>();
    }
    x[i][j] = k;
}

class ErrorRateTable {
public:
    ErrorRateTable(uint n_qubits);
    ErrorRateTable(uint n_qubits, fp_t def_p);  // Initialize all tables to same value.
    ErrorRateTable(uint n_qubits, fp_t def_1q_p, def_2q_p);  // Initializes single
                                                             // and two-qubit tables
                                                             // respectively.
    //      Entries:    Opname      Affected pair          P(II), P(IX), ..., P(ZZ)
    typedef std::tuple<std::string, std::pair<uint, uint>, std::array<fp_t, 16>> 
        corr_t;

    map2d<std::string, uint, fp_t> op1q;
    map2d<std::string, std::pair<uint, uint>, fp_t> op2q;

    // Other errors: not set by default.
    map2d<std::string, std::pair<uint, uint>, fp_t> op2q_leakage;
    map2d<std::string, std::pair<uint, uint>, fp_t> op2q_crosstalk;
    map2d<std::string, std::pair<uint, uint>, std::vector<corr_t>> op2q_correlated;
private:
    void initialize(void);

    const uint n_qubits;
    fp_t def_1q = 0.0;
    fp_t def_2q = 0.0;
};

class TimeTable {  // Units are in nanoseconds.
public:
    TimeTable(uint n_qubits);
    // Initializes the tables to the default times provided.
    TimeTable(uint n_qubits, fp_t def_op, fp_t def_t1, fp_t def_t2);
    TimeTable(uint n_qubits, fp_t def_1q, fp_t def_2q, fp_t def_t1, fp_t def_t2);

    map2d<std::string, uint, fp_t> op1q;
    map2d<std::string, std::pair<uint, uint>, fp_t> op2q;
    std::map<uint, fp_t> t1;
    std::map<uint, fp_t> t2;
private:
    void initialize(void);

    const uint n_qubits;
    fp_t def_1q = 30;
    fp_t def_2q = 40;
    fp_t def_t1 = 15000;
    fp_t def_t2 = 7500;

    const fp_t def_ro = 620;
};

} // qontra

#endif // TABLES_h
