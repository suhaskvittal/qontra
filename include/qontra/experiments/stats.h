/*
 *  author: Suhas Vittal
 *  date:   15 December 2023
 * */

#ifndef QONTRA_EXPERIMENTS_STATS_h
#define QONTRA_EXPERIMENTS_STATS_h

#include <mpi.h>
#include <stdint.h>

namespace qontra {

template <class T> MPI_Datatype get_mpi_type(void);

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
    static statistic_t from_range(MPI_Op, PTR, PTR, size_t);

    void    reduce(void);

    // Python style slicing:
    statistic_t slice(size_t from, size_t to);

    void        fill(T);
    size_t      size(void) const;
    void        resize(size_t);

    std::vector<T>  vec(void) const;

    // Accesses to data in this class is via statistic_t[i] -- gets the i-th element.
    //
    // [] returns a reference, at() returns a value.
    inline T& operator[](size_t index);
    inline T at(size_t index=0) const;

    // Arithmetic operations below (some of which are overloads).
    //
    // Note that these are not well-defined for all types.
    void    negate(void);
    void    root(void);
    // Scalar arithmetic operations (operation is applied on all elements).
    statistic_t&    operator+=(T);
    statistic_t&    operator-=(T);
    statistic_t&    operator*=(T);
    statistic_t&    operator/=(T);
    // Vector arithmetic operations (operation is applied elementwise).
    statistic_t&    operator+=(statistic_t&);
    statistic_t&    operator-=(statistic_t&);
    statistic_t&    operator*=(statistic_t&);
    statistic_t&    operator/=(statistic_t&);
    // Comparison and update operations.
    void                        scalar_replace_if_better_extrema(T);
    template <class PTR> void   vector_replace_if_better_extrema(PTR array);

    // Some basic statistic functions:
    //      For functions that rely on another statistic_t object, if the sizes do not line up,
    //      then the output is meaningless.
    statistic_t<fp_t>   get_mean(uint64_t n);
    statistic_t<fp_t>   get_std(statistic_t<fp_t> mean, uint64_t n);

    template <class LAMBDA> void                apply_scalar_elementwise(LAMBDA);
    template <class LAMBDA, class PTR> void     apply_vector_elementwise(PTR other, LAMBDA);

    const MPI_Datatype  mpi_type;
    const MPI_Op        mpi_op;
private:
    T compare_and_return_extrema(T, T);

    std::vector<T> data;
};

}   // qontra

#include "stats.inl"

#endif  // QONTRA_EXPERIMENTS_STATS_h
