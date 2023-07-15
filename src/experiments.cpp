/*
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#include "experiments.h"

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

void
generate_syndromes(const stim::Circuit& circuit, 
                    uint64_t shots, 
                    callback_t callbacks)
{
    uint64_t __shots;

    int world_size, world_rank;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        __shots = shots / world_size;
        if (world_rank == 0)    __shots += shots % world_size;
    } else {
        world_size = 1;
        world_rank = 0;
        __shots = shots;
    }

    std::mt19937_64 rng(G_BASE_SEED + world_rank);

    while (__shots > 0) {
        const uint64_t shots_this_batch = 
            __shots < G_SHOTS_PER_BATCH ? __shots : G_SHOTS_PER_BATCH;
        auto samples = stim::detector_samples(
                            circuit, shots_this_batch, false, true, rng);
        samples = samples.transposed();
        for (uint64_t t = 0; t < shots_this_batch; t++) {
            stim::simd_bits_range_ref row = samples[t];
            callbacks.prologue(row);
        }
        __shots -= shots_this_batch;
    }
}

void
build_syndrome_trace(std::string output_folder,
        const stim::Circuit& circuit, 
        uint64_t shots) 
{
    int world_size = 1, world_rank = 0;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    uint width = circuit.count_detectors() + circuit.count_observables();
    stim::simd_bit_table syndromes(G_SHOTS_PER_BATCH, width);

    syndromes.clear();

    stim::simd_bits ref(width);

    uint64_t ctr = 0;
    uint file_offset = world_rank;
    cb_t1 t_cb = [&] (stim::simd_bits_range_ref row) {
        syndromes[++ctr] |= row;
        if (ctr == G_SHOTS_PER_BATCH) {
            ctr = 0;
            std::string filename = output_folder + "/batch_"
                                + std::to_string(file_offset)
                                + ".dets";
            auto output_trace = syndromes.transposed();
            FILE* fout = fopen(filename.c_str(), "w");
            stim::write_table_data(fout,
                                    G_SHOTS_PER_BATCH,
                                    width,
                                    ref,
                                    output_trace,
                                    stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                    'D',
                                    'L',
                                    circuit.count_detectors());
            fclose(fout);
            syndromes.clear();
            file_offset += world_rank;
        }
    };
    callback_t srcbs;
    srcbs.prologue = t_cb;
    generate_syndromes(circuit, shots, srcbs);
    // Write out any remaining syndromes as well.
    if (ctr) {
        std::string filename = output_folder + "/batch_"
                            + std::to_string(file_offset)
                            + ".dets";
        auto output_trace = syndromes.transposed();
        FILE* fout = fopen(filename.c_str(), "w");
        stim::write_table_data(fout,
                                G_SHOTS_PER_BATCH,
                                width,
                                ref,
                                output_trace,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                'L',
                                circuit.count_detectors());
        fclose(fout);
    }
}

uint64_t
read_syndrome_trace(std::string folder, 
                        const stim::Circuit& circuit,
                        callback_t callbacks) 
{
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();

    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }

    uint file_offset = world_rank;
    std::string batch_file = folder + "/batch_" 
                                + std::to_string(file_offset) + ".dets";
    uint64_t local_shots = 0;
    while (access(batch_file.c_str(), F_OK) >= 0) {
        stim::simd_bit_table samples(G_SHOTS_PER_BATCH, det+obs);
        FILE* fin = fopen(batch_file.c_str(), "r");
        uint64_t true_shots = stim::read_file_data_into_shot_table(
                                fin,
                                G_SHOTS_PER_BATCH,
                                det,
                                stim::SampleFormat::SAMPLE_FORMAT_DETS,
                                'D',
                                samples,
                                true,
                                0,
                                det,
                                obs);
        for (uint64_t s = 0; s < true_shots; s++) {
            stim::simd_bits_range_ref row = samples[s];
            callbacks.prologue(row);
        }
        file_offset += world_size;
        batch_file = folder + "/batch_" + std::to_string(file_offset) + ".dets";
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
memory_experiment(decoder::Decoder* dec, experiments::memory_params_t params) {
    const stim::Circuit circuit = dec->get_circuit();
    // Stats that will be placed in memory_result_t
    uint64_t    logical_errors, __logical_errors = 0;
    uint64_t    hw_sum, __hw_sum = 0;
    uint64_t    hw_sqr_sum, __hw_sqr_sum = 0;
    uint64_t    hw_max, __hw_max = 0;
    fp_t        t_sum, __t_sum = 0;
    fp_t        t_sqr_sum, __t_sqr_sum = 0;
    fp_t        t_max, __t_max = 0;

    const uint sample_width = 
        circuit.count_detectors() + circuit.count_observables();
    cb_t1 dec_cb = [&] (stim::simd_bits_range_ref& row)
    {
        params.callbacks.prologue(row);
        uint obs_bits = 0;
        for (uint i = 0; i < circuit.count_observables(); i++) {
            obs_bits += row[i+circuit.count_detectors()];
        }
        const uint hw = row.popcnt() - obs_bits;
        __hw_sum += hw;
        __hw_sqr_sum += SQR(hw);
        __hw_max = hw > __hw_max ? hw : __hw_max;
        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
            return;
        }
        syndrome_t syndrome(row);
        auto res = dec->decode_error(syndrome); 
        __logical_errors += res.is_error;
        __t_sum += res.exec_time;
        __t_sqr_sum += SQR(res.exec_time);
        __t_max = res.exec_time > __t_max ? res.exec_time : __t_max;
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

    // Collect results using MPI if enabled.
    if (G_USE_MPI) {
        MPI_Reduce(&__logical_errors, &logical_errors, 1, MPI_UNSIGNED_LONG,
                    MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&__hw_sum, &hw_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                    0, MPI_COMM_WORLD);
        MPI_Reduce(&__hw_sqr_sum, &hw_sqr_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                    0, MPI_COMM_WORLD);
        MPI_Reduce(&__hw_max, &hw_max, 1, MPI_UNSIGNED_LONG, MPI_MAX,
                    0, MPI_COMM_WORLD);
        MPI_Reduce(&__t_sum, &t_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&__t_sqr_sum, &t_sqr_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 
                    0, MPI_COMM_WORLD);
        MPI_Reduce(&__t_max, &t_max, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    } else {
        logical_errors = __logical_errors;
        hw_sum = __hw_sum;
        hw_sqr_sum = __hw_sqr_sum;
        hw_max = __hw_max;
        t_sum = __t_sum;
        t_sqr_sum = __t_sqr_sum;
        t_max = __t_max;
    }
    fp_t ler = ((fp_t)logical_errors) / ((fp_t)shots);
    fp_t hw_mean = MEAN(hw_sum, shots);
    fp_t hw_std = STD(hw_mean, hw_sqr_sum, shots);
    fp_t t_mean = MEAN(t_sum, shots);
    fp_t t_std = STD(t_mean, t_sqr_sum, shots);

    memory_result_t res = {
        ler,
        hw_mean,
        hw_std,
        hw_max,
        t_mean,
        t_std,
        t_max
    };
    return res;
}

}   // qontra
