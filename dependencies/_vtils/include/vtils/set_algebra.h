/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef VTILS_SET_ALGEBRA_h
#define VTILS_SET_ALGEBRA_h

namespace vtils {

template <class SET_TYPE> SET_TYPE immd_set_union(SET_TYPE, SET_TYPE);
template <class SET_TYPE> SET_TYPE immd_set_diff(SET_TYPE, SET_TYPE);
template <class SET_TYPE> SET_TYPE immd_set_intersect(SET_TYPE, SET_TYPE);
template <class SET_TYPE> SET_TYPE immd_set_symdiff(SET_TYPE, SET_TYPE);

template <class SET_TYPE> SET_TYPE& upd_set_union(SET_TYPE&, SET_TYPE);
template <class SET_TYPE> SET_TYPE& upd_set_diff(SET_TYPE&, SET_TYPE);
template <class SET_TYPE> SET_TYPE& upd_set_intersect(SET_TYPE&, SET_TYPE);
template <class SET_TYPE> SET_TYPE& upd_set_symdiff(SET_TYPE&, SET_TYPE);

}   // vtils

#include "inl/set_algebra.inl"

#endif  // VTILS_SET_ALGEBRA_h
