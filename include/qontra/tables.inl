/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <limits>

namespace qontra {

inline ErrorTable&
ErrorTable::operator=(const ErrorTable& other) {
    op1q = other.op1q;
    op2q = other.op2q;
    idling = other.idling;
    m1w0 = other.m1w0;
    m0w1 = other.m0w1;
    op2q_leakage_injection = other.op2q_leakage_injection;
    op2q_leakage_transport = other.op2q_leakage_transport;
    return *this;
}

inline TimeTable&
TimeTable::operator=(const TimeTable& other) {
    op1q = other.op1q;
    op2q = other.op2q;
    t1 = other.t1;
    t2 = other.t2;
    return *this;
}

inline ErrorTable
operator*(ErrorTable t, fp_t x) {
    t *= x; 
    return t;
}

inline ErrorTable
operator*(fp_t x, ErrorTable t) {
    return t*x;
}

inline TimeTable
operator*(TimeTable t, fp_t x) {
    t *= x;
    return t;
}

inline TimeTable
operator*(fp_t x, TimeTable t) {
    return t*x;
}

namespace tables {

inline ErrorAndTiming
operator*(ErrorAndTiming et, fp_t x) {
    et *= x; return et;
}

inline ErrorAndTiming
operator*(fp_t x, ErrorAndTiming et) {
    return et * x;
}

}   // tables

}   // qontra
