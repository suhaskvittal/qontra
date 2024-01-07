/*
 *  author: Suhas Vittal
 *  date:   2 September 2023
 * */

#include "mldh/latency_sim.h"

namespace qontra {

using namespace experiments;

namespace mldh {

latency_sim_stats_t
simulate_on_latency_data(std::string data_folder, const latency_sim_params_t& params) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    Manager mgr(data_folder);

    // Normalized latency = true_syndromes / total_syndromes.
    uint64_t total_syndromes = 0, __total_syndromes = 0;
    uint64_t true_syndromes = 0, __true_syndromes = 0;
    
    // We need to track, for every syndrome, how many rounds
    // have elapsed since we received it, and how much
    // time will we take (a running total).
    uint64_t rounds = 1;
    uint64_t total_time = 0;

    uint64_t t;
    while (mgr.get_next(&t)) {
        uint64_t tt = t + (rounds-1)*params.round_latency;
        // tt is just t, but relative to the original
        // syndrome. If tt > total, then our new
        // running total is tt, as decoding this
        // syndrome is now the bottleneck.
        total_time = tt > total_time ? tt : total_time;

        if (total_time <= rounds*params.round_latency) {
            // We have settled the backlog.
            __true_syndromes++;
            // Reset total and rounds.
            total_time = 0;
            rounds = 1;
        } else {
            rounds++;
        }
        __total_syndromes++;
    }
    // Compute stats.
    latency_sim_stats_t stats;
    if (G_USE_MPI) {
        MPI_Allreduce(&__total_syndromes, &total_syndromes, 1,
                    MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&__true_syndromes, &true_syndromes, 1,
                    MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    } else {
        total_syndromes = __total_syndromes;
        true_syndromes = __true_syndromes;
    }
    stats.normalized_latency = ((fp_t) total_syndromes) / ((fp_t) true_syndromes);

    return stats;
}

}   // mldh
}   // qontra
