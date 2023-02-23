/*
 *  author: Suhas Vittal
 *  date:   22 February 2023
 * */

#include "benchmark/statbench/numerical_dist.h"

namespace qrc {
namespace benchmark {
namespace statbench {

//#define NUMERICAL_IS_NORMAL

NumericalDistribution::NumericalDistribution(const NumericalModelData& model_data)
    :ErrorDistribution(model_data.sample_mean, model_data.sample_variance, true),
    model(model_data.model),
    log_sample_mean(log(model_data.sample_mean)),
    log_sample_variance(log(model_data.sample_variance)),
    sample_size(model_data.sample_size)
{}

fp_t
NumericalDistribution::pmf(fp_t hw, bool logscale) {
    fp_t hwdiv2 = hw * 0.5;

    fp_t log_prob = 1.0;
    if (model == NumericalModel::normal) {
        log_prob = -0.5*(log(2*M_PI) + log(variance)) - 0.5*(hwdiv2-mean)*(hwdiv2-mean)/variance;
    } else if (model == NumericalModel::lognormal) {
        log_prob = -log(hwdiv2) - 0.5*(log(2*M_PI) + log(variance)) 
                    - 0.5 * (log(hwdiv2)-mean)*(log(hwdiv2)-mean)/variance;
    }
    
    if (logscale) {
        return log_prob;
    } else {
        return pow(M_E, log_prob);
    }
}

fp_t
NumericalDistribution::cdf(fp_t hw, bool logscale) {
    fp_t hwdiv2 = hw * 0.5;
    // Convert to z-score
    fp_t sx = sqrt(variance / sample_size);
    fp_t z = (hwdiv2 - mean) / sx;
    fp_t log_prob = 0.5 + log(2.0 - erfc(z/sqrt(2)));

    if (logscale) {
        return log_prob;
    } else {
        return pow(M_E, log_prob);
    }
}

NumericalDistribution
build_numerical_distribution(const std::map<uint, uint64_t>& hamming_weight_dist, NumericalModel model) {
    fp_t sample_mean = 0.0;
    fp_t sample_size = 0.0;

    fp_t m2 = 0.0;
    for (auto pair : hamming_weight_dist) {
        sample_size += pair.second;
        fp_t n = (fp_t) pair.second;
        fp_t x = ceil( ((fp_t)pair.first) * 0.5);

        fp_t delta = n * (x - sample_mean);
        sample_mean += delta / sample_size;
        m2 += delta*(x - sample_mean);
    }

    fp_t sample_variance = m2 / (sample_size - 1.0);
    NumericalModelData data = {sample_mean, sample_variance, sample_size};
    NumericalDistribution dist(data);
    return dist;
}

}   // statbench
}   // benchmark
}   // qrc
