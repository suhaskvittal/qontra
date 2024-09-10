/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#ifndef VTILS_UTILITY_h
#define VTILS_UTILITY_h

#include <initializer_list>

namespace vtils {

template <class VECTOR_LIKE, class CONTAINER>
void push_back_range(VECTOR_LIKE&, CONTAINER, int begin_offset=0, int end_offset=0);

template <class VECTOR_LIKE, class ITER>
void push_back_range(VECTOR_LIKE&, ITER begin, ITER end);

template <class VECTOR_LIKE, class T>
void push_back_all(VECTOR_LIKE&, std::initializer_list<T>);

template <class SET_LIKE, class CONTAINER>
void insert_range(SET_LIKE&, CONTAINER, int begin_offset=0, int end_offset=0);

template <class SET_LIKE, class ITER>
void insert_range(SET_LIKE&, ITER begin, ITER end);

template <class SET_LIKE, class T>
void insert_all(SET_LIKE&, std::initializer_list<T>);

}   // vtils

#include "inl/utility.inl"

#endif  // VTILS_UTILITY_h
