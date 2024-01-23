/*
 *  author: Suhas Vittal
 *  date    23 January 2024
 * */

#include "qontra/experiments.h"

#include <iostream>
#include <map>
#include <set>
#include <string>

#include <math.h>
#include <string.h>

namespace qontra {

template <class T> inline MPI_Datatype get_mpi_type() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype get_mpi_type<uint64_t>() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype get_mpi_type<int64_t>() { return MPI_LONG; }
template <> inline MPI_Datatype get_mpi_type<double>() { return MPI_DOUBLE; }

template <class T> template <class PTR> inline statistic_t<T>
statistic_t<T>::from_range(MPI_Op op, PTR start, PTR end, size_t max_size) {
    statistic_t s(op, max_size);
    size_t true_size = 0;
    while (start != end) {
        s[true_size++] = *start;
        start++;
    }
    s.resize(true_size);
    return s;
}

template <class T> void
statistic_t<T>::reduce() {
    if (!G_USE_MPI) return;

    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size == 1) return;    // Nothing to be done.

    std::vector<T> buf(size());
    MPI_Allreduce(&data[0], &buf[0], size(), mpi_type, mpi_op, MPI_COMM_WORLD);
    data = std::move(buf);
}

template <class T> inline statistic_t<T>
statistic_t<T>::slice(size_t from, size_t to) {
    T* start = &data[from];
    T* end = start + (to-from);
    return from_range(mpi_op, start, end, to-from);
}

template <class T> inline void              statistic_t<T>::fill(T x) { memset(&data[0], x, size()); }
template <class T> inline size_t            statistic_t<T>::size() const { return data.size(); }
template <class T> inline void              statistic_t<T>::resize(size_t s) { data.resize(s); }
template <class T> inline std::vector<T>    statistic_t<T>::vec() const { return data; }

template <class T> inline T&    statistic_t<T>::operator[](size_t index) { return data[index]; }
template <class T> inline T     statistic_t<T>::at(size_t index) const { return data.at(index); }

template <class T> inline void
statistic_t<T>::negate(void) {
    apply_scalar_elementwise([] (T y) { return -y; }); 
}

template <class T> inline void
statistic_t<T>::root(void) {
    apply_scalar_elementwise([] (T y) { return sqrt(y); });
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator+=(T x) {
    apply_scalar_elementwise([&] (T y) { return y + x; }); 
    return *this; 
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator-=(T x)  { 
    apply_scalar_elementwise([&] (T y) { return y + x; }); 
    return *this;
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator*=(T x) {
    apply_scalar_elementwise([&] (T y) { return y * x; }); 
    return *this;
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator/=(T x) {
    apply_scalar_elementwise([&] (T y) { return y / x; });
    return *this; 
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator+=(statistic_t<T>& other) { 
    apply_vector_elementwise(other, [] (T x, T y) { return x + y; });
    return *this;
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator-=(statistic_t<T>& other) { 
    apply_vector_elementwise(other, [] (T x, T y) { return x - y; });
    return *this;
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator*=(statistic_t<T>& other) { 
    apply_vector_elementwise(other, [] (T x, T y) { return x * y; });
    return *this;
}

template <class T> inline statistic_t<T>&
statistic_t<T>::operator/=(statistic_t<T>& other) { 
    apply_vector_elementwise(other, [] (T x, T y) { return x / y; });
    return *this;
}

template <class T> inline void
statistic_t<T>::scalar_replace_if_better_extrema(T x) {
    apply_scalar_elementwise([&] (T y) { return compare_and_return_extrema(x, y); });
}

template <class T> template <class PTR> inline void
statistic_t<T>::vector_replace_if_better_extrema(PTR array) {
    apply_vector_elementwise(array, [&] (T x, T y) { return compare_and_return_extrema(x, y); });
}

template <class T> inline statistic_t<fp_t>
statistic_t<T>::get_mean(uint64_t n) {
    statistic_t<fp_t> mean(*this);
    mean /= static_cast<fp_t>(n);
    return mean;
}

template <class T> inline statistic_t<fp_t>
statistic_t<T>::get_std(statistic_t<fp_t> mean, uint64_t n) {
    statistic_t<fp_t> stddev(*this);
    if (mean.size() != size()) return stddev;
    // First compute variance.
    stddev /= static_cast<fp_t>(n);
    stddev -= mean;
    // Now take the elementwise square root.
    stddev.root();
    return stddev;
}

template <class T> template <class LAMBDA> inline void
statistic_t<T>::apply_scalar_elementwise(LAMBDA f) {
    for (size_t i = 0; i < size(); i++) {
        data[i] = f(data[i]);
    }
}

template <class T> template <class LAMBDA, class PTR> inline void
statistic_t<T>::apply_vector_elementwise(PTR other, LAMBDA f) {
    for (size_t i = 0; i < size(); i++) {
        data[i] = f(data[i], static_cast<T>(other[i]));
    }
}

template <class T> inline T
statistic_t<T>::compare_and_return_extrema(T x, T y) {
    if (mpi_op == MPI_MAX) {
        return x > y ? x : y;
    } else {
        return x < y ? x : y;
    }
}

}   // qontra
