/*
 *  author: Suhas Vittal
 *  date:   22 February 2023
 * */

#ifndef BENCHMARK_STATBENCH_ANALYTICAL_DIST_h
#define BENCHMARK_STATBENCH_ANALYTICAL_DIST_h

#include <stim.h>

#include "benchmark/statbench.h"
#include "benchmark/statbench/numerical_dist.h"
#include "decoding_graph.h"
#include "defs.h"

#include <math.h>

namespace qrc {
namespace benchmark {
namespace statbench {

enum class AnalyticalModel {binomial, poisson};

struct AnalyticalModelData {
    fp_t mean;
    fp_t variance;
    fp_t n_events;
    fp_t error_rate;
    
    AnalyticalModel model;
};

class AnalyticalDistribution : public ErrorDistribution {
public:
    AnalyticalDistribution(const AnalyticalModelData&);

    fp_t pmf(fp_t, bool) override;
    fp_t cdf(fp_t, bool) override;

    AnalyticalModel model;
private:
    // Parameters for distributions.
    fp_t n_events;
    fp_t error_rate;
    fp_t log_error_rate;
    fp_t log_error_rate_complement;

    fp_t log_mean;
};

AnalyticalDistribution
build_analytical_distribution(const stim::Circuit&, AnalyticalModel);
AnalyticalDistribution
bootstrap_from_normal_approximation(const stim::Circuit&, AnalyticalModel, const NumericalDistribution&);

}   // statbench
}   // benchmark
}   // qrc

#endif  // BENCHMARK_STATBENCH_ANALYTICAL_DIST_h
