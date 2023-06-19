/*
 *  author: Suhas Vittal
 *  date:   12 June 2023
 * */

#include "experiments.h"

#define SQR(x)              (x)*(x)
#define MEAN(s, n)          ((fp_t)(s))/((fp_t)(n))
#define STD(m, ss, n)       ( ((fp_t)(ss))/((fp_t)(n)) - SQR(m) )

namespace qontra {

namespace experiments {

bool        G_USE_MPI  = true;
uint64_t    G_SHOTS_PER_BATCH = 100'000;
uint64_t    G_BASE_SEED = 0;
bool        G_FILTER_OUT_SYNDROMES = true;
uint64_t    G_FILTERING_HAMMING_WEIGHT = 2;

}   // experiments

using namespace experiments;

void
generate_syndromes(const stim::Circuit& circuit, uint64_t shots, cb_t1 cb) {
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
            cb(row);
        }
        __shots -= shots_this_batch;
    }
}

memory_result_t
memory_experiment(decoder::Decoder* dec, uint64_t shots, cb_t1 cb1, cb_t2 cb2) {
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
        cb1(row);
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
        cb2(res);
    };

    generate_syndromes(circuit, shots, dec_cb);

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
