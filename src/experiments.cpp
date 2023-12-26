/*
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#include "experiments.h"
#include "stats.h"

#include <fcntl.h>
#include <mpi.h>
#include <stdio.h>
#include <strings.h>
#include <unistd.h>

namespace qontra {

namespace experiments {

callback_t  DEFAULT_CALLBACKS;

bool        G_USE_MPI  = true;
uint64_t    G_SHOTS_PER_BATCH = 100'000;
uint64_t    G_BASE_SEED = 0;
bool        G_FILTER_OUT_SYNDROMES = true;
uint64_t    G_FILTERING_HAMMING_WEIGHT = 2;

}   // experiments

using namespace experiments;

template <class T> T sqr(T x) { return x*x; }

inline std::string
get_batch_filename(uint batchno) {
    std::string batch_file = std::string("/batch_")
                                + std::to_string(batchno)
                                + std::string(".dets");
    return batch_file;
}

inline void
configure_optimal_batch_size() {
    uint64_t cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);  // In bytes.
    G_SHOTS_PER_BATCH = cache_line_size << 3;   // Need to convert to bits.
}

void
generate_syndromes(const stim::Circuit& circuit, uint64_t shots, callback_t callbacks) {
    uint64_t local_shots = shots;

    int world_size = 1, world_rank = 0;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        local_shots = shots / world_size;
        if (world_rank == 0) local_shots += shots % world_size;
    }

    std::mt19937_64 rng(G_BASE_SEED + world_rank);

    stim::DetectorErrorModel dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(
                                        circuit,
                                        true,   // decompose errors
                                        true,   // fold loops
                                        true,   // allow gauge detectors
                                        0.0,
                                        false,  // ignore decomposition failures
                                        true    // block decomposition from introducing remnant edges
                                    );
    stim::DemSampler<SIMD_WIDTH> sampler(dem, std::move(rng), G_SHOTS_PER_BATCH);
    while (local_shots > 0) {
        const uint64_t shots_this_batch = local_shots < G_SHOTS_PER_BATCH ? local_shots : G_SHOTS_PER_BATCH;
        sampler.resample(false);

        stim::simd_bit_table<SIMD_WIDTH> syndrome_table = sampler.det_buffer.transposed(),
                                            obs_table = sampler.obs_buffer.transposed();
        for (size_t t = 0; t < shots_this_batch; t++) {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = syndrome_table[t],
                                                    obs = obs_table[t];
            callbacks.prologue({syndrome, obs});
        }
        local_shots -= shots_this_batch;
    }
}

void
build_syndrome_trace(std::string output_folder, const stim::Circuit& circuit, uint64_t shots) {
    int world_size = 1, world_rank = 0;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    const uint w_det = circuit.count_detectors();
    const uint w_obs = circuit.count_observables();

    stim::simd_bit_table<SIMD_WIDTH> syndromes(G_SHOTS_PER_BATCH, w_det);
    stim::simd_bit_table<SIMD_WIDTH> observables(G_SHOTS_PER_BATCH, w_obs);

    syndromes.clear();

    stim::simd_bits<SIMD_WIDTH> ref(w_det+w_obs);
    ref.clear();

    uint64_t ctr = 0;
    uint file_offset = world_rank;

    // Create a local function for writing syndromes and observables to a file.
    auto write_batch = [&] (void) {
        ctr = 0;
        // Merge syndrome and observable tables into one.
        stim::simd_bit_table<SIMD_WIDTH> output_trace(std::move(syndromes));
        output_trace = output_trace.transposed();
        output_trace.resize(w_det+w_obs, G_SHOTS_PER_BATCH);

        stim::simd_bit_table<SIMD_WIDTH> obs_tr(std::move(observables));
        obs_tr = obs_tr.transposed();
        for (uint i = 0; i < w_obs; i++) output_trace[w_det + i] |= obs_tr[i];

        std::string filename = output_folder + "/" + get_batch_filename(file_offset);

        FILE* fout = fopen(filename.c_str(), "w");
        stim::write_table_data(fout,
                                G_SHOTS_PER_BATCH,
                                w_det+w_obs,
                                ref,
                                output_trace,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                'L',
                                w_det);
        fclose(fout);
        syndromes.clear();
        observables.clear();
        file_offset += world_rank;
    };

    cb_t1 t_cb = [&] (shot_payload_t payload) {
        // Record payload to simd_bit_tables.
        syndromes[ctr] |= payload.syndrome;
        observables[ctr] |= payload.observables;
        ctr++;

        if (ctr == G_SHOTS_PER_BATCH) {
            write_batch();
        }
    };
    callback_t srcbs;
    srcbs.prologue = t_cb;
    generate_syndromes(circuit, shots, srcbs);
    // Write out any remaining syndromes as well.
    if (ctr) {
        write_batch();
    }
}

uint64_t
read_syndrome_trace(std::string folder, const stim::Circuit& circuit, callback_t callbacks) {
    const uint w_det = circuit.count_detectors();
    const uint w_obs = circuit.count_observables();

    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    uint file_offset = world_rank;
    uint64_t local_shots = 0;

    std::string batch_file = folder + "/" + get_batch_filename(file_offset);
    while (faccessat(AT_FDCWD, batch_file.c_str(), F_OK, 0) >= 0) {
        // We will temporarily write the data column-wise to input_trace, and then extract
        // the syndrome and observable tables from there.
        stim::simd_bit_table<SIMD_WIDTH> input_trace(w_det+w_obs, G_SHOTS_PER_BATCH);

        FILE* fin = fopen(batch_file.c_str(), "r");
        uint64_t true_shots = stim::read_file_data_into_shot_table(
                                fin,
                                G_SHOTS_PER_BATCH,
                                w_det+w_obs,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                input_trace,
                                false);
        fclose(fin);
        input_trace.destructive_resize(w_det+w_obs, true_shots);
        // Split the data into two tables.
        stim::simd_bit_table<SIMD_WIDTH> syndromes(std::move(input_trace));
        stim::simd_bit_table<SIMD_WIDTH> observables(w_obs, true_shots);
        // First populate observables.
        for (uint i = 0; i < w_obs; i++) observables[i] |= syndromes[w_det+i];
        observables = observables.transposed();
        // Now, we will create syndromes by destructive resizing + transposing.
        syndromes.destructive_resize(w_det, true_shots);
        syndromes = syndromes.transposed();
        // Finally, execute the callbacks.
        for (uint64_t s = 0; s < true_shots; s++) {
            stim::simd_bits_range_ref<SIMD_WIDTH> syndrome = syndromes[s],
                                                    obs = observables[s];
            callbacks.prologue({syndrome, obs});
        }
        file_offset += world_size;
        batch_file = folder + "/" + get_batch_filename(file_offset);
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

memory_result_t
memory_experiment(Decoder* dec, experiments::memory_params_t params) {
    const stim::Circuit circuit = dec->get_circuit();
    const uint w_det = circuit.count_detectors();
    const uint w_obs = circuit.count_observables();

    statistic_t<uint64_t> logical_errors(MPI_SUM, w_obs+1);
    statistic_t<uint64_t> hw_sum(MPI_SUM), hw_sqr_sum(MPI_SUM), hw_max(MPI_MAX);
    statistic_t<fp_t> t_sum(MPI_SUM), t_sqr_sum(MPI_SUM), t_max(MPI_MAX);

    const uint sample_width = w_det + w_obs;
    cb_t1 dec_cb = [&] (shot_payload_t payload)
    {
        stim::simd_bits<SIMD_WIDTH> syndrome(payload.syndrome),
                                    obs(payload.observables);

        params.callbacks.prologue({syndrome, obs});
        const uint hw = syndrome.popcnt();
        // Update HW statistics and skip the trial if the HW is too small
        // and filtering is enabled.
        hw_sum += hw;
        hw_sqr_sum += sqr(hw);
        hw_max.scalar_replace_if_better_extrema(hw);

        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
            return;
        }
        // Decode syndrome
        auto res = dec->decode_error(syndrome); 
        logical_errors[0] += (bool) (payload.observables != res.corr);
        t_sum += res.exec_time;
        t_sqr_sum += sqr(res.exec_time);
        t_max.scalar_replace_if_better_extrema(res.exec_time);

        for (uint i = 0; i < w_obs; i++) {
            logical_errors[i+1] += (bool)(res.corr[i] != obs[i]);
        }

        params.callbacks.epilogue(res);
    };

    callback_t srcbs;
    srcbs.prologue = dec_cb;

    uint64_t shots = params.shots;
    if (shots == 0) {
        shots = read_syndrome_trace(params.trace_folder, circuit, srcbs);
    } else {
        generate_syndromes(circuit, shots, srcbs);
    }

    // Collect results across all processors.
    logical_errors.reduce();

    hw_sum.reduce();
    hw_sqr_sum.reduce();
    hw_max.reduce();

    t_sum.reduce();
    t_sqr_sum.reduce();
    t_max.reduce();
    // Compute means and variances/std. deviations.
    statistic_t<fp_t> logical_error_rate = logical_errors.get_mean(shots);
    statistic_t<fp_t> hw_mean = hw_sum.get_mean(shots);
    statistic_t<fp_t> t_mean = t_sum.get_mean(shots);

    statistic_t<fp_t> hw_std = hw_sqr_sum.get_std(hw_mean, shots);
    statistic_t<fp_t> t_std = t_sqr_sum.get_std(t_mean, shots);

    std::vector<fp_t> logical_error_rate_by_obs = logical_error_rate.slice(1, w_obs+1).vec();

    memory_result_t res = {
        logical_error_rate.at(),
        hw_mean.at(),
        hw_std.at(),
        hw_max.at(),
        t_mean.at(),
        t_std.at(),
        t_max.at(),
        // Additional statistics:
        logical_error_rate_by_obs
    };
    return res;
}

}   // qontra
