/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <random>

#include <fcntl.h>
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

namespace qontra {

inline void
configure_optimal_batch_size() {
#ifndef L1D_CACHE_LINE_SIZE
#if defined(linux)
    uint64_t cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);  // In bytes.
#else
    uint64_t cache_line_size = 64;
#endif
#else
    uint64_t cache_line_size = L1D_CACHE_LINE_SIZE;
#endif
    G_SHOTS_PER_BATCH = cache_line_size << 3;   // Need to convert to bits.
}

inline std::string
get_batch_filename(size_t batchno) {
    std::string batch_file = std::string("/batch_")
                                + std::to_string(batchno)
                                + std::string(".dets");
    return batch_file;
}

template <class CALLBACK> void
generate_syndromes(const DetailedStimCircuit& circuit, uint64_t shots, CALLBACK cb) {
    uint64_t local_shots = shots;

    int world_size = 1, world_rank = 0;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        local_shots = shots / world_size + (world_rank == 0) * (shots % world_size);
    }

    std::mt19937_64 rng(G_BASE_SEED + world_rank);
    while (local_shots > 0) {
        const uint64_t shots_this_batch = local_shots < G_SHOTS_PER_BATCH ? local_shots : G_SHOTS_PER_BATCH;
        auto sample = stim::sample_batch_detection_events<SIMD_WIDTH>(circuit, shots_this_batch, rng);
        stim::simd_bit_table<SIMD_WIDTH> syndrome_table = sample.first.transposed(),
                                            obs_table = sample.second.transposed();

        for (size_t t = 0; t < shots_this_batch; t++) {
            cb({syndrome_table[t], obs_table[t]});
        }
        local_shots -= shots_this_batch;
    }
}

template <class CALLBACK> uint64_t
read_syndrome_trace(std::string input_folder, const DetailedStimCircuit& circuit, CALLBACK cb) {
    const size_t n_det = circuit.count_detectors();
    const size_t n_obs = circuit.count_observables();

    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    int file_offset = world_rank;
    uint64_t local_shots = 0;

    std::string batch_file = input_folder + "/" + get_batch_filename(file_offset);
    while (faccessat(AT_FDCWD, batch_file.c_str(), F_OK, 0) >= 0) {
        // We will temporarily write the data column-wise to input_trace, and then extract
        // the syndrome and observable tables from there.
        stim::simd_bit_table<SIMD_WIDTH> input_trace(n_det+n_obs, G_SHOTS_PER_BATCH);

        FILE* fin = fopen(batch_file.c_str(), "r");
        uint64_t true_shots = stim::read_file_data_into_shot_table(
                                fin,
                                G_SHOTS_PER_BATCH,
                                n_det+n_obs,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                input_trace,
                                false);
        fclose(fin);
        input_trace.destructive_resize(n_det+n_obs, true_shots);
        // Split the data into two tables.
        stim::simd_bit_table<SIMD_WIDTH> syndromes(std::move(input_trace));
        stim::simd_bit_table<SIMD_WIDTH> observables(n_obs, true_shots);
        // First populate observables.
        for (size_t i = 0; i < n_obs; i++) {
            observables[i] |= syndromes[n_det+i];
        }
        observables = observables.transposed();
        // Now, we will create syndromes by destructive resizing + transposing.
        syndromes.destructive_resize(n_det, true_shots);
        syndromes = syndromes.transposed();
        // Finally, execute the callbacks.
        for (uint64_t s = 0; s < true_shots; s++) {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = syndromes[s],
                                                    obs = observables[s];
            cb({syndromes[s], observables[s]});
        }
        file_offset += world_size;
        batch_file = input_folder + "/" + get_batch_filename(file_offset);
        local_shots += true_shots;
    }
    uint64_t shots;
    if (G_USE_MPI) {
        MPI_Allreduce(&local_shots, &shots, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                        MPI_COMM_WORLD);
    } else {
        shots = local_shots;
    }
    return shots;
}

}   // qontra
