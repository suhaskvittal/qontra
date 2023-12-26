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
        data(size)
    {
        fill(default_value);
    }

    statistic_t(const statistic_t& other)
        :mpi_type(other.mpi_type),
        mpi_op(other.mpi_op),
        data(other.data)
    {}

    template <class U>
    statistic_t(const statistic_t<U>& other)
        :mpi_type(get_mpi_type<T>()),
        mpi_op(other.mpi_op),
        data(other.size())
    {
        // Cast the other statistic_t's data pointers to the type for this object.
        for (size_t i = 0; i < size(); i++) {
            data[i] = static_cast<T>(other.at(i));
        }
    }

    statistic_t(statistic_t&& other)
        :mpi_type(other.mpi_type),
        mpi_op(other.mpi_op),
        data(std::move(other.data))
    {}

    template <class PTR>
    static statistic_t from_range(MPI_Op op, PTR start, PTR end, size_t max_size) {
        statistic_t s(op, max_size);
        size_t true_size = 0;
        while (start != end) {
            s[true_size++] = *start;
            start++;
        }
        s.resize(true_size);
        return s;
    }

    void reduce() {
        int world_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        if (world_size == 1) return;    // Nothing to be done.

        std::vector<T> buf(size());
        MPI_Allreduce(&data[0], &buf[0], size(), mpi_type, mpi_op, MPI_COMM_WORLD);
        data = std::move(buf);
    }

    // Python style slicing:
    inline statistic_t slice(size_t from, size_t to) {
        T* start = &data[from];
        T* end = start + (to-from);
        return from_range(mpi_op, start, end, to-from);
    }

    inline void fill(T x) { memset(&data[0], x, size()); }
    inline size_t size() const { return data.size(); }
    inline void resize(size_t s) { data.resize(s); }

    inline std::vector<T> vec(void) { return data; }

    // Accesses to data in this class is via statistic_t[i] -- gets the i-th element.
    //
    // [] returns a reference, at() returns a value.
    inline T& operator[](int index) { return data[index]; }
    inline T at(int index=0) const { return data[index]; }

    // Arithmetic operations below (some of which are overloads).
    //
    // Note that these are not well-defined for all types.
    inline void negate(void) { apply_scalar_elementwise([] (T y) { return -y; }); }
    inline void root(void) { apply_scalar_elementwise([] (T y) { return sqrt(y); }); }

    inline statistic_t& operator+=(T x) 
        { apply_scalar_elementwise([&] (T y) { return y + x; }); return *this; }
    inline statistic_t& operator-=(T x) 
        { apply_scalar_elementwise([&] (T y) { return y + x; }); return *this; }
    inline statistic_t& operator*=(T x) 
        { apply_scalar_elementwise([&] (T y) { return y * x; }); return *this; }
    inline statistic_t& operator/=(T x) 
        { apply_scalar_elementwise([&] (T y) { return y / x; }); return *this; }

    inline statistic_t& operator+=(statistic_t& other) { 
        apply_vector_elementwise(other, [] (T x, T y) { return x + y; });
        return *this;
    }

    inline statistic_t& operator-=(statistic_t& other) { 
        apply_vector_elementwise(other, [] (T x, T y) { return x - y; });
        return *this;
    }

    inline statistic_t& operator*=(statistic_t& other) { 
        apply_vector_elementwise(other, [] (T x, T y) { return x * y; });
        return *this;
    }

    inline statistic_t& operator/=(statistic_t& other) { 
        apply_vector_elementwise(other, [] (T x, T y) { return x / y; });
        return *this;
    }

    inline void scalar_replace_if_better_extrema(T x) {
        apply_scalar_elementwise([&] (T y) { return compare_and_return_extrema(x, y); });
    }

    template <class PTR>
    inline void vector_replace_if_better_extrema(PTR array) {
        apply_vector_elementwise(array, [&] (T x, T y) { return compare_and_return_extrema(x, y); });
    }

    // Some basic statistic functions:
    //      For functions that rely on another statistic_t object, if the sizes do not line up,
    //      then the output is meaningless.

    inline statistic_t<fp_t> get_mean(uint64_t n) {
        statistic_t<fp_t> mean(*this);
        mean /= static_cast<fp_t>(n);
        return mean;
    }

    statistic_t<fp_t> get_std(statistic_t<fp_t> mean, uint64_t n) {
        statistic_t<fp_t> stddev(*this);
        if (mean.size() != size()) return stddev;
        // First compute variance.
        stddev /= static_cast<fp_t>(n);
        stddev -= mean;
        // Now take the elementwise square root.
        stddev.root();
        return stddev;
    }

    template <class LAMBDA> inline void
    apply_scalar_elementwise(LAMBDA f) {
        for (size_t i = 0; i < size(); i++) {
            data[i] = f(data[i]);
        }
    }

    template <class LAMBDA, class PTR> inline void
    apply_vector_elementwise(PTR other, LAMBDA f) {
        for (size_t i = 0; i < size(); i++) {
            data[i] = f(data[i], static_cast<T>(other[i]));
        }
    }

    const MPI_Datatype  mpi_type;
    const MPI_Op        mpi_op;
private:
    inline T compare_and_return_extrema(T x, T y) {
        if (mpi_op == MPI_MAX) {
            return x > y ? x : y;
        } else {
            return x < y ? x : y;
        }
    }

    std::vector<T> data;
};

}   // qontra

#endif  // QONTRA_STATS_h
