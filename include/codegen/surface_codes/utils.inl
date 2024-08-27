/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

namespace cgen {

template <int R, int S, class T, int M> bool
Check<R,S,T,M>::add_qubit(sptr<Qubit<R,S>> q) {
    if (is_full()) return false;
    qubits[size++] = q;
    _cast()->update_qubit(q);
    return true;
}

template <int R, int S> inline void
Star<R,S>::tie(Star<R,S>* other, uint8_t i, uint8_t j) {
    const sptr<Qubit<R,S>>& q = this->at(i);
    // First, delete other[j].
#ifdef SC_DEBUG
    // Check that q and other[j] is otherwise empty (z_ptr == 1 and x_ptr == 0).
    for (sptr<Qubit<R,S>> x : {q,other[j]}) {
        if (x->z_ptr != 1 || x->x_ptr != 0) {
            std::cerr << "[ Star::tie ] found z_ptr = " << x->z_ptr
                << " and x_ptr = " << x->x_ptr << " during tie.\n";
        }
    }
#endif
    other->qubits[j] = q;
    other->update_qubit(q);
}

template <int R, int S> inline bool
Star<R,S>::update_qubit(sptr<Qubit<R,S>> q) {
    if (q->x_ptr == 2) return false;
    q->x_checks[q->x_ptr++] = this;
    return true;
}

}   // cgen
