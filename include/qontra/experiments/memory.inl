/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include "qontra/decoder.h"
#include "qontra/experiments.h"
#include "qontra/experiments/stats.h"
#include "qontra/ext/qes.h"

namespace qontra {

inline qontra::DetailedStimCircuit
make_circuit(std::string qes_file, fp_t p) {
    qes::Program<> program = qes::from_file(qes_file);
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;
    et = et * (1000*p);

    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing);
}

inline memory_result_t
memory_experiment(Decoder* dec, memory_config_t config) {
    return memory_experiment(dec, config, [] (shot_payload_t x) {}, [] (Decoder::result_t x) {});
}

template <class PROLOGUE, class EPILOGUE> memory_result_t
memory_experiment(Decoder* dec, memory_config_t config, PROLOGUE p_cb, EPILOGUE e_cb) {
    const DetailedStimCircuit circuit = dec->get_circuit();
    const size_t n_obs = circuit.count_observables();

    statistic_t<uint64_t> logical_errors(MPI_SUM, n_obs+1);
    statistic_t<uint64_t> hw_sum(MPI_SUM), hw_sqr_sum(MPI_SUM), hw_max(MPI_MAX);
    statistic_t<fp_t> t_sum(MPI_SUM), t_sqr_sum(MPI_SUM), t_max(MPI_MAX);

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

        if (experiments::G_FILTER_OUT_SYNDROMES
                && hw <= experiments::G_FILTERING_HAMMING_WEIGHT) 
        {
            return;
        }
        // Decode syndrome
        auto res = dec->decode_error(syndrome); 
        logical_errors[0] += (bool) (payload.observables != res.corr);
        t_sum += res.exec_time;
        t_sqr_sum += sqr(res.exec_time);
        t_max.scalar_replace_if_better_extrema(res.exec_time);

        for (uint i = 0; i < n_obs; i++) {
            logical_errors[i+1] += (bool)(res.corr[i] != obs[i]);
        }
        e_cb(res);
    };

    uint64_t shots = config.shots;
    if (shots == 0) {
        shots = read_syndrome_trace(config.trace_folder, circuit, dec_cb);
    } else {
        generate_syndromes(circuit, shots, dec_cb);
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
