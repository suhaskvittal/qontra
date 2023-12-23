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

#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>

namespace qontra {

struct statistic_t;
//
// Declare operator overloads first.
//
statistic_t operator-(statistic_t); 

statistic_t operator+(statistic_t, uint64_t);
statistic_t operator-(statistic_t, uint64_t);
statistic_t operator*(statistic_t, uint64_t);
statistic_t operator/(statistic_t, uint64_t);

statistic_t operator+(statistic_t, int64_t);
statistic_t operator-(statistic_t, int64_t);
statistic_t operator*(statistic_t, int64_t);
statistic_t operator/(statistic_t, int64_t);

statistic_t operator+(statistic_t, double);
statistic_t operator-(statistic_t, double);
statistic_t operator*(statistic_t, double);
statistic_t operator/(statistic_t, double);
//
// Now declare rest of the struct.
//
struct statistic_t {
    statistic_t(void)
        :name(""),
        mpi_type(MPI_UNSIGNED_LONG),
        mpi_op(MPI_SUM),
        u64(0),
        i64(0),
        fp64(0)
    {}

    statistic_t(std::string name, MPI_Datatype mpi_type, MPI_Op mpi_op)
        :name(name),
        mpi_type(mpi_type),
        mpi_op(mpi_op),
        u64(0),
        i64(0),
        fp64(0)
    {}

    statistic_t(const statistic_t& other)
        :name(other.name),
        mpi_type(other.mpi_type),
        mpi_op(other.mpi_op),
        u64(other.u64),
        i64(other.i64),
        fp64(other.fp64)
    {}

    std::string name;

    MPI_Datatype    mpi_type;
    MPI_Op          mpi_op;
    // Any other type can be downcast from these.
    uint64_t    u64;
    int64_t     i64;
    double      fp64;

    std::string str(void) {
        if (mpi_type == MPI_UNSIGNED_LONG)  return std::to_string(u64);
        else if (mpi_type == MPI_LONG)      return std::to_string(i64);
        else                                return std::to_string(fp64);
    }

#define STAT_OO_UPDATE(op, type) statistic_t& operator op##= (type x) { *this = *this op x; return *this; }
    
    STAT_OO_UPDATE(+, uint64_t)
    STAT_OO_UPDATE(-, uint64_t)
    STAT_OO_UPDATE(*, uint64_t)
    STAT_OO_UPDATE(/, uint64_t)

    STAT_OO_UPDATE(+, int64_t)
    STAT_OO_UPDATE(-, int64_t)
    STAT_OO_UPDATE(*, int64_t)
    STAT_OO_UPDATE(/, int64_t)

    STAT_OO_UPDATE(+, double)
    STAT_OO_UPDATE(-, double)
    STAT_OO_UPDATE(*, double)
    STAT_OO_UPDATE(/, double)
};

class StatFile {
public:
    StatFile() 
        :backing_store()
    {}

    void decl(std::string name, MPI_Datatype mpi_type, MPI_Op mpi_op) {
        backing_store[name] = statistic_t(name, mpi_type, mpi_op); 
        variables.insert(name);
    }

    void remove(std::string name) { backing_store.erase(name); }
    statistic_t& operator[](std::string x) { return backing_store[x]; }
    void clear(void) { backing_store.clear(); }

    std::set<std::string> get_variables(void) { return variables; }

    StatFile mpi_acc(void);
private:
    std::set<std::string> variables;
    std::map<std::string, statistic_t>  backing_store;
};

}   // qontra

#endif  // QONTRA_STATS_h
