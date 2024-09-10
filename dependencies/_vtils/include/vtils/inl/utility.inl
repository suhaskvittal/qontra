/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

namespace vtils {

template <class VECTOR_LIKE, class CONTAINER> inline void
push_back_range(VECTOR_LIKE& x, CONTAINER y, int boff, int eoff) {
    x.insert(x.end(), std::next(y.begin(), boff), std::next(y.end(), eoff));
}

template <class VECTOR_LIKE, class ITER> inline void
push_back_range(VECTOR_LIKE& x, ITER begin, ITER end) {
    x.insert(x.end(), begin, end);
}

template <class VECTOR_LIKE, class T> inline void
push_back_all(VECTOR_LIKE& x, std::initializer_list<T> args) {
    for (T y : args) x.push_back(y);
}

template <class SET_LIKE, class CONTAINER> inline void
insert_range(SET_LIKE& x, CONTAINER y, int boff, int eoff) {
    x.insert(std::next(y.begin(), boff), std::next(y.end(), eoff));
}

template <class SET_LIKE, class ITER> inline void
insert_range(SET_LIKE& x, ITER begin, ITER end) {
    x.insert(begin, end);
}

template <class SET_LIKE, class T> inline void
insert_all(SET_LIKE& x, std::initializer_list<T> args) {
    for (T y : args) x.insert(y);
}

}   // vtils
