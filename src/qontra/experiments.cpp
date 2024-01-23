/*
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#include "qontra/experiments.h"

#include <mpi.h>
#include <stdio.h>

namespace qontra {

bool        G_USE_MPI  = true;
uint64_t    G_SHOTS_PER_BATCH = 100'000;
uint64_t    G_BASE_SEED = 0;
bool        G_FILTER_OUT_SYNDROMES = true;
uint64_t    G_FILTERING_HAMMING_WEIGHT = 2;

void
build_syndrome_trace(std::string output_folder, const stim::Circuit& circuit, uint64_t shots) {
    int world_size = 1, world_rank = 0;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    const size_t n_det = circuit.count_detectors();
    const size_t n_obs = circuit.count_observables();

    stim::simd_bit_table<SIMD_WIDTH> syndromes(G_SHOTS_PER_BATCH, n_det);
    stim::simd_bit_table<SIMD_WIDTH> observables(G_SHOTS_PER_BATCH, n_obs);

    syndromes.clear();

    stim::simd_bits<SIMD_WIDTH> ref(n_det+n_obs);
    ref.clear();

    uint64_t ctr = 0;
    int file_offset = world_rank;

    // Create a local function for writing syndromes and observables to a file.
    auto write_batch = [&] (void) {
        ctr = 0;
        // Merge syndrome and observable tables into one.
        stim::simd_bit_table<SIMD_WIDTH> output_trace(std::move(syndromes));
        output_trace = output_trace.transposed();
        output_trace.resize(n_det+n_obs, G_SHOTS_PER_BATCH);

        stim::simd_bit_table<SIMD_WIDTH> obs_tr(std::move(observables));
        obs_tr = obs_tr.transposed();
        for (size_t i = 0; i < n_obs; i++) {
            output_trace[n_det + i] |= obs_tr[i];
        }

        std::string filename = output_folder + "/" + get_batch_filename(file_offset);

        FILE* fout = fopen(filename.c_str(), "w");
        stim::write_table_data(fout,
                                G_SHOTS_PER_BATCH,
                                n_det+n_obs,
                                ref,
                                output_trace,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                'L',
                                n_det);
        fclose(fout);
        syndromes.clear();
        observables.clear();
        file_offset += world_rank;
    };

    generate_syndromes(circuit, shots,
            [&] (shot_payload_t payload) {
                // Record payload to simd_bit_tables.
                syndromes[ctr] = std::move(payload.syndrome);
                observables[ctr] = std::move(payload.observables);
                ctr++;

                if (ctr == G_SHOTS_PER_BATCH) {
                    write_batch();
                }
            });
    // Write out any remaining syndromes as well.
    if (ctr) {
        write_batch();
    }
}

}   // qontra
