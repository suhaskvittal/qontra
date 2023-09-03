/*
 *  author: Suhas Vittal
 *  date:   2 September 2023
 * */

#ifndef MLDH_LATENCY_SIM_h
#define MLDH_LATENCY_SIM_h

#include "defs.h"
#include "experiments.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include <math.h>
#include <mpi.h>

namespace qontra {
namespace mldh {

struct latency_sim_params_t {
    uint64_t round_latency = 1000;

    bool simulate_parallel_decoder = true;
};

struct latency_sim_stats_t {
    fp_t normalized_latency;
};

// This latency simulatiion uses the timing
// data from a real decoder to simulate its
// expected slowdown in a real program.

latency_sim_stats_t
simulate_on_latency_data(std::string data_folder, const latency_sim_params_t&);

// This class is a simple Manager for retrieving
// times from the data folder.
class Manager {
public:
    Manager(std::string data_folder)
        :data_folder(data_folder),
        eod(false),
        fin(),
        fno(0),
        mpi_world_rank(0),
        mpi_world_size(1)
    {
        if (experiments::G_USE_MPI) {
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);

            fin = std::ifstream(get_kth_file(fno));
        }
    }

    // Returns false if there is no more data.
    bool get_next(uint64_t* t_p) {
        if (eod)    return false;

        fin.read(reinterpret_cast<char*>(t_p), sizeof(uint64_t));
        if (fin.eof()) {
            fno++;
            // Open next file.
            std::filesystem::path input_path(get_kth_file(fno));
            if (!std::filesystem::exists(input_path)) {
                eod = true;
            }
            fin = std::ifstream(input_path);
            fin.seekg(mpi_world_rank*sizeof(uint64_t));
        } else {
            fin.seekg((mpi_world_size-1)*sizeof(uint64_t), std::ios::cur);
        }
        return true;
    }

    int mpi_world_rank;
    int mpi_world_size;
private:
    std::string get_kth_file(uint64_t k) {
        return data_folder + std::string("/proc_") + std::to_string(k) + std::string(".bin");
    }

    const std::string data_folder;
    bool eod;

    std::ifstream fin;
    uint64_t fno;
};

}   // mldh
}   // qontra

#endif  // MLDH_LATENCY_SIM_h
