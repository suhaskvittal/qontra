/*
 *  author: Suhas Vittal
 *  date:   3 February 2023
 * */

#ifndef BENCHMARK_STATBENCH_h
#define BENCHMARK_STATBENCH_h

#include "defs.h"

#include <functional>
#include <map>

namespace qrc {

class Decoder;
/*
 *  Decoder pointers from a dgf_t should be
 *  allocated on the heap.
 * */
typedef std::function<Decoder*(fp_t)> dgf_t;

namespace benchmark {

struct StatisticalResult {
    fp_t n_logical_errors = 0;
    fp_t mean_execution_time = 0;
    fp_t max_execution_time = 0;
    uint64_t true_shots = 0;
    fp_t statistical_shots = 0;
    std::map<uint, fp_t> hamming_weight_dist;
};

class ErrorDistribution {
public:
    ErrorDistribution(fp_t mean, fp_t variance, bool is_cdf_supported)
        :mean(mean), variance(variance), is_cdf_supported(is_cdf_supported) {}

    virtual fp_t pmf(fp_t hw, bool logscale) =0;
    virtual fp_t cdf(fp_t hw, bool logscale) =0;

    fp_t get_mean(void) const { return mean; }
    fp_t get_variance(void) const { return variance; }

    bool cdf_supported(void) { return is_cdf_supported; }
protected:
    fp_t mean;
    fp_t variance;

    bool is_cdf_supported;
};

}   // benchmark
}   // qrc

#endif // BENCHMARK_STATBENCH
