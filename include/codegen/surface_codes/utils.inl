/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */
#ifdef SC_DEBUG
#include <iostream>
#endif

namespace cgen {

template <int R, int S, class T, int M> bool
Check<R,S,T,M>::add_qubit(sptr<Qubit<R,S>> q) {
    if (is_full()) return false;
    qubits[size++] = q;
    _cast()->update_qubit(q);
    return true;
}

template <int R, int S> inline void
Star<R,S>::tie(Star<R,S>* other) {
    sptr<Qubit<R,S>> q = std::make_shared<Qubit<R,S>>();
    this->add_qubit(q);
    other->add_qubit(q);
}

template <int R, int S> inline bool
Star<R,S>::update_qubit(sptr<Qubit<R,S>> q) {
    if (q->x_ptr == 2) return false;
    q->x_checks[q->x_ptr++] = this;
    return true;
}

}   // cgen
