/*
 *  author: Suhas Vittal
 *  date:   22 February 2023
 * */

#ifndef BENCHMARK_STATBENCH_NUMERICAL_DIST_h
#define BENCHMARK_STATBENCH_NUMERICAL_DIST_h

#include "benchmark/statbench.h"
#include "defs.h"

#include <math.h>

namespace qrc {
namespace benchmark {
namespace statbench {

enum class NumericalModel {normal, lognormal};

struct NumericalModelData {
    fp_t sample_mean;
    fp_t sample_variance;
    fp_t sample_size;

    NumericalModel model;
};

class NumericalDistribution : public ErrorDistribution {
public:
    NumericalDistribution(const NumericalModelData&);

    fp_t pmf(fp_t, bool) override;
    fp_t cdf(fp_t, bool) override;

    NumericalModel model;
private:
    fp_t log_sample_mean;
    fp_t log_sample_variance;
    fp_t sample_size;
};

NumericalDistribution
build_numerical_distribution(const std::map<uint, uint64_t>& hamming_weight_dist, NumericalModel);

}   // statbench
}   // benchmark
}   // qrc

#endif  // BENCHMARK_STATBENCH_NUMERICAL_DIST_h
