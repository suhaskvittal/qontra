/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include "qontra/experiments.h"
#include "qontra/experiments/stats.h"
#include "qontra/ext/qes.h"

#ifdef DECODER_PERF
#include <vtils/timer.h>
#endif

namespace qontra {

inline DetailedStimCircuit
make_default_circuit(std::string qes_file, fp_t p, bool fix_timing_error_as_p, std::string model) {
    qes::Program<> program = qes::from_file(qes_file);
    return make_default_circuit(program, p, fix_timing_error_as_p, model);
}

DetailedStimCircuit
make_default_circuit(
        const qes::Program<>& program, fp_t p, bool fix_timing_error_as_p, std::string model) 
{
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;
    et = et * (1000*p);

    if (model == "cap") {
        et.e_m1w0 = 0.0;
        et.e_m0w1 = 0.0;
        et.e_g1q = 0.0;
        et.e_g2q = 0.0;
        et.e_idle = 0.0;
    } else if (model == "pheno") {
        et.e_g1q = 0.0;
        et.e_g2q = 0.0;
        et.e_idle = 0.0;
    }

    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing, fix_timing_error_as_p ? p : -1.0);
}

inline memory_result_t
run_memory_with_generated_syndromes(Decoder* dec, memory_config_t config) {
    return run_memory_with_generated_syndromes(
            dec, config, [] (shot_payload_t x) {}, [] (Decoder::result_t x) {});
}

template <class PROLOGUE, class EPILOGUE> memory_result_t
run_memory_with_generated_syndromes(Decoder* dec, memory_config_t config, PROLOGUE p_cb, EPILOGUE e_cb) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }

    const DetailedStimCircuit circuit = dec->get_circuit();
    const size_t n_obs = circuit.count_observables();

    statistic_t<uint64_t> logical_errors(MPI_SUM, n_obs+1);
    statistic_t<uint64_t> hw_sum(MPI_SUM), hw_sqr_sum(MPI_SUM), hw_max(MPI_MAX);
    statistic_t<fp_t> t_sum(MPI_SUM), t_sqr_sum(MPI_SUM), t_max(MPI_MAX);

    statistic_t<uint64_t> shots(MPI_SUM);

    auto dec_cb = [&] (shot_payload_t payload)
    {
        stim::simd_bits<SIMD_WIDTH> syndrome(payload.syndrome),
                                    obs(payload.observables);
        p_cb({syndrome, obs});
        const size_t hw = syndrome.popcnt();
        // Update HW statistics and skip the trial if the HW is too small
        // and filtering is enabled.
        hw_sum += hw;
        hw_sqr_sum += sqr(hw);
        hw_max.scalar_replace_if_better_extrema(hw);
        shots += 1;

        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT && obs.popcnt() == 0) {
            return;
        }
#ifdef MEMORY_DEBUG
        std::cout << "--------------------------------------" << std::endl;
#ifdef DECODER_PERF
        vtils::Timer timer;
        fp_t t;
        timer.clk_start();
#endif
#endif
        // Decode syndrome
        auto res = dec->decode_error(syndrome); 
#ifdef MEMORY_DEBUG
#ifdef DECODER_PERF
        t = timer.clk_end();
        std::cout << "[ memory ] decoder took " << t*1e-9 << "s to decode HW " << hw << std::endl;
#endif
        if (payload.observables != res.corr) {
            std::cout << "is logical error!" << std::endl;
            std::cout << "\texpected: ";
            for (size_t i = 0; i < n_obs; i++) std::cout << payload.observables[i]+0;
            std::cout << std::endl << "\treceived: ";
            for (size_t i = 0; i < n_obs; i++) std::cout << res.corr[i]+0;
            std::cout << std::endl;
        }
#endif
        logical_errors[0] += (bool) (payload.observables != res.corr);
        t_sum += res.exec_time;
        t_sqr_sum += sqr(res.exec_time);
        t_max.scalar_replace_if_better_extrema(res.exec_time);

        for (size_t i = 0; i < n_obs; i++) {
            logical_errors[i+1] += (bool)(res.corr[i] != obs[i]);
        }
        e_cb(res);
    };

    const uint64_t orig_base_seed = G_BASE_SEED;

    uint64_t n_errors = 0;
    while (n_errors < config.errors_until_stop) {
        generate_syndromes(circuit, world_size*G_SHOTS_PER_BATCH, dec_cb);
        G_BASE_SEED += world_size;
        if (G_USE_MPI) {
            MPI_Allreduce(&logical_errors[0], &n_errors, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        } else {
            n_errors = logical_errors[0];
        }
    }
    G_BASE_SEED = orig_base_seed;

    // Collect results across all processors.
    logical_errors.reduce();

    hw_sum.reduce();
    hw_sqr_sum.reduce();
    hw_max.reduce();

    t_sum.reduce();
    t_sqr_sum.reduce();
    t_max.reduce();

    shots.reduce();
    uint64_t _shots = shots.at();
    if (world_rank == 0) {
        std::cout << "found " << logical_errors.at() << " errors within " << _shots << " shots" << std::endl;
    }
    // Compute means and variances/std. deviations.
    statistic_t<fp_t> logical_error_rate = logical_errors.get_mean(_shots);
    statistic_t<fp_t> hw_mean = hw_sum.get_mean(_shots);
    statistic_t<fp_t> t_mean = t_sum.get_mean(_shots);

    statistic_t<fp_t> hw_std = hw_sqr_sum.get_std(hw_mean, _shots);
    statistic_t<fp_t> t_std = t_sqr_sum.get_std(t_mean, _shots);

    std::vector<fp_t> logical_error_rate_by_obs = logical_error_rate.slice(1, n_obs+1).vec();

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
