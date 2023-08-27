/*
 *  author: Suhas Vittal
 *  date:   27 August 2023
 * */

#ifndef QONTRA_EXT_CYTNX
#define QONTRA_EXT_CYTNX

#include <Scalar.hpp>
#include <Tensor.hpp>

namespace qontra {
    
template <typename T> cytnx::Tensor
operator*(cytnx::Tensor t, T x) {
    return t * cytnx::Scalar(x, Type.ComplexDouble);
}

template <typename T> cytnx::Tensor
operator*(T x, cytnx::Tensor t) {
    return t * x;
}

}   // qontra

#endif  // QONTRA_CYTNX_EXT
