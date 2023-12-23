/*
 *  author: Suhas Vittal
 *  date:   15 December 2023
 * */

#include "stats.h"

namespace qontra {

StatFile
StatFile::mpi_acc() {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size == 1) return *this;

    StatFile acc_stats;

    for (auto& pair : backing_store) {
        std::string name = pair.first;
        statistic_t& st = pair.second;
        // Accumulate result.
        void* data_p;
        if (st.mpi_type == MPI_UNSIGNED_LONG) {
            data_p = &st.u64;
        } else if (st.mpi_type == MPI_LONG) {
            data_p = &st.i64;
        } else if (st.mpi_type == MPI_DOUBLE) {
            data_p = &st.fp64;
        } else {
            std::cerr << "StatFile::mpi_acc : unsupported data type (MPI type = " 
                    << st.mpi_type << ")\n"
                    << "Only types supported are: MPI_UNSIGNED_LONG, "
                    << "MPI_DOUBLE, and MPI_LONG." << std::endl;
            exit(1);
        }

        void* buf = malloc(sizeof(uint64_t));  // All data types are 8 bytes.
        MPI_Allreduce(data_p, buf, 1, st.mpi_type, MPI_SUM, MPI_COMM_WORLD);
        // Create accumulated statistic.
        acc_stats.decl(name, st.mpi_type, st.mpi_op);
        acc_stats[name].u64 = *(uint64_t*)buf;
        acc_stats[name].i64 = *(int64_t*)buf;
        acc_stats[name].fp64 = *(double*)buf;
        free(buf);
    }
    return acc_stats;
}

// Operator overloads:
#define __STAT_OO_U64(op) statistic_t operator op (statistic_t st, uint64_t x)\
                            {\
                                st.u64 op## = x;\
                                st.i64 op## = x;\
                                st.fp64 op## = x;\
                                return st;\
                            }

#define __STAT_OO_I64(op) statistic_t operator op (statistic_t st, int64_t x)\
                            {\
                                if (st.mpi_type == MPI_UNSIGNED_LONG) st.mpi_type = MPI_LONG;\
                                st.u64 op## = x;\
                                st.i64 op## = x;\
                                st.fp64 op## = x;\
                                return st;\
                            }

#define __STAT_OO_FP64(op) statistic_t operator op (statistic_t st, double x)\
                            {\
                                st.mpi_type = MPI_DOUBLE;\
                                st.fp64 op## = x;\
                                return st;\
                            }

statistic_t operator-(statistic_t st) {
    st.i64 = -st.i64;
    st.fp64 = -st.fp64;
    return st;
}

__STAT_OO_U64(+)
__STAT_OO_U64(-)
__STAT_OO_U64(*)
__STAT_OO_U64(/)

__STAT_OO_I64(+)
__STAT_OO_I64(-)
__STAT_OO_I64(*)
__STAT_OO_I64(/)

__STAT_OO_FP64(+)
__STAT_OO_FP64(-)
__STAT_OO_FP64(*)
__STAT_OO_FP64(/)


}   // qontra
