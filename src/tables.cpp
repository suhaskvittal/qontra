/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#include "tables.h"

namespace contra {

ErrorRateTable::ErrorRateTable(uint n_qubits)
    :ErrorRateTable(n_qubits, 0, 0)
{ }

ErrorRateTable::ErrorRateTable(uint n_qubits, fp_t def_p)
    :ErrorRateTable(n_qubits, def_p, def_p)
{ }

ErrorRateTable::ErrorRateTable(uint n_qubits, fp_t def_1q_p, fp_t def_2q_p)
    :n_qubits(n_qubits),
    def_1q(def_1q_p),
    def_2q(def_2q_p)
{
    initialize();
}

void
ErrorRateTable::initialize() {
    for (uint i = 0; i < n_qubits; i++) {
        put(op1q, "H", i, def_1q);
        put(op1q, "X", i, def_1q);
        put(op1q, "Z", i, def_1q);
        put(op1q, "M", i, def_1q);
        put(op1q, "R", i, def_1q);

        for (uint j = 0; j < n_qubits; j++) {
            if (i == j) continue;
            auto i_j = std::make_pair(i, j);
            put(op2q, "CX", i_j, def_2q);
            put(op2q_leakage, "CX", i_j, 0);
            put(op2q_crosstalk, "CX", i_j, 0);
            put(op2q_correlated, "CX", i_j, std::vector<corr_t>());
        }
    }
}

TimeTable::TimeTable(uint n_qubits)
    :TimeTable(n_qubits, 30, 40, 15000, 7500)
{ }

TimeTable::TimeTable(uint n_qubits, fp_t def_op, fp_t def_t1, fp_t def_t2)
    :TimeTable(n_qubits, def_op, def_op, def_t1, def_t2)
{ }

TimeTable::TimeTable(uint n_qubits, fp_t def_1q, fp_t def_2q, fp_t def_t1, fp_t def_t2)
    :TimeTable(n_qubits, def_1q, def_2q, def_t1, def_t2)
{
    initialize();
}

void
TimeTable::initialize() {
    for (uint i = 0; i < n_qubits; i++) {
        put(op1q, "H", i, def_1q);
        put(op1q, "X", i, def_1q);
        put(op1q, "Z", i, def_1q);
        put(op1q, "M", i, def_ro);
        put(op1q, "R", i, def_1q);
        t1[i] = def_t1;
        t2[i] = def_t2;

        for (uint j = 0; j < n_qubits; j++) {
            if (i == j) continue;
            auto i_j = std::make_pair(i, j);
            put(op2q, "CX", i_j, def_2q);
        }
    }
}

}   // contra
