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
    op2q_crosstalk = other.op2q_crosstalk;
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

namespace tables {

inline ErrorAndTiming&
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
