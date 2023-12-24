/*
 *  author: Suhas Vittal
 *  date:   15 December 2023
 * */

#ifndef QONTRA_STATS_h
#define QONTRA_STATS_h

#include <iostream>
#include <map>
#include <set>
#include <string>

#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <string.h>

namespace qontra {

template <class T> inline MPI_Datatype get_mpi_type() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype get_mpi_type<uint64_t>() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype get_mpi_type<int64_t>() { return MPI_LONG; }
template <> inline MPI_Datatype get_mpi_type<double>() { return MPI_DOUBLE; }

template <class T>
class statistic_t {
public:
    statistic_t(MPI_Op op, int size=1, T default_value=0)
        :mpi_type(get_mpi_type<T>()),
        mpi_op(op),
        data_p(std::make_unique<T[]>(size)),
        size(size)
    {
        memset(data_p.get(), default_value, size); // Clear out memory.
    }

    statistic_t(const statistic_t& other)
        :mpi_type(other.mpi_type),
        mpi_op(other.mpi_op),
        data_p(std::make_unique<T[]>(other.size)),
        size(other.size)
    {
        memcpy(data_p.get(), other.data_p.get(), size);
    }

    template <class U>
    statistic_t(const statistic_t<U>& other)
        :mpi_type(get_mpi_type<T>()),
        mpi_op(other.mpi_op),
        data_p(std::make_unique<T[]>(other.size)),
        size(other.size)
    {
        // Cast the other statistic_t's data pointers to the type for this object.
        for (int i = 0; i < size; i++) {
            data_p[i] = static_cast<T>(other(i));
        }
    }

    statistic_t(statistic_t&& other)
        :mpi_type(other.mpi_type),
        mpi_op(other.mpi_op),
        data_p(std::move(other.data_p)),
        size(other.size)
    {}

    void reduce() {
        int world_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        if (world_size == 1) return;    // Nothing to be done.

        uptr<T[]> buf = std::make_unique<T[]>(size);
        MPI_Allreduce(data_p.get(), buf.get(), size, mpi_type, mpi_op, MPI_COMM_WORLD);
        data_p = std::move(buf);
    }

    // Accesses to data in this class is via statistic_t[i] -- gets the i-th element.
    //
    // [] returns a reference, () returns a value.
    T& operator[](int index) {
        return data_p[index];
    }

    T operator()(int index=0) const {
        return data_p[index];
    }

    // Arithmetic operations below (some of which are overloads).
    //
    // Note that these are not well-defined for all types.
    void negate(void) { apply_scalar_elementwise([] (T y) { return -y; }); }
    void root(void) { apply_scalar_elementwise([] (T y) { return sqrt(y); }); }

    statistic_t& operator+=(T x) { apply_scalar_elementwise([&] (T y) { return y + x; }); return *this; }
    statistic_t& operator-=(T x) { apply_scalar_elementwise([&] (T y) { return y + x; }); return *this; }
    statistic_t& operator*=(T x) { apply_scalar_elementwise([&] (T y) { return y * x; }); return *this; }
    statistic_t& operator/=(T x) { apply_scalar_elementwise([&] (T y) { return y / x; }); return *this; }

    statistic_t& operator+=(statistic_t& other) { 
        apply_vector_elementwise(other.data_p.get(), [] (T x, T y) { return x + y; });
        return *this;
    }

    statistic_t& operator-=(statistic_t& other) { 
        apply_vector_elementwise(other.data_p.get(), [] (T x, T y) { return x - y; });
        return *this;
    }

    statistic_t& operator*=(statistic_t& other) { 
        apply_vector_elementwise(other.data_p.get(), [] (T x, T y) { return x * y; });
        return *this;
    }

    statistic_t& operator/=(statistic_t& other) { 
        apply_vector_elementwise(other.data_p.get(), [] (T x, T y) { return x / y; });
        return *this;
    }

    void scalar_replace_if_better_extrema(T x) {
        apply_scalar_elementwise([&] (T y) { return compare_and_return_extrema(x, y); });
    }

    template <class PTR>
    void vector_replace_if_better_extrema(PTR array) {
        apply_vector_elementwise(array, [&] (T x, T y) { return compare_and_return_extrema(x, y); });
    }

    // Some basic statistic functions:
    //      For functions that rely on another statistic_t object, if the sizes do not line up,
    //      then the output is meaningless.

    statistic_t<fp_t> get_mean(uint64_t n) {
        statistic_t<fp_t> mean(*this);
        mean /= static_cast<fp_t>(n);
        return mean;
    }

    statistic_t<fp_t> get_std(statistic_t<fp_t> mean, uint64_t n) {
        statistic_t<fp_t> stddev(*this);
        if (mean.size != size) return stddev;
        // First compute variance.
        stddev /= static_cast<fp_t>(n);
        stddev -= mean;
        // Now take the elementwise square root.
        stddev.root();
        return stddev;
    }

    template <class LAMBDA> void
    apply_scalar_elementwise(LAMBDA f) {
        for (int i = 0; i < size; i++) {
            data_p[i] = f(data_p[i]);
        }
    }

    template <class LAMBDA, class PTR> void
    apply_vector_elementwise(PTR other, LAMBDA f) {
        for (int i = 0; i < size; i++) {
            data_p[i] = f(data_p[i], static_cast<T>(other[i]));
        }
    }

    const MPI_Datatype  mpi_type;
    const MPI_Op        mpi_op;
    const int size;
private:
    T compare_and_return_extrema(T x, T y) {
        if (mpi_op == MPI_MAX) {
            return x > y ? x : y;
        } else {
            return x < y ? x : y;
        }
    }

    uptr<T[]> data_p;
};

}   // qontra

#endif  // QONTRA_STATS_h
