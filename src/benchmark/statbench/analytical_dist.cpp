/*
 *  author: Suhas Vittal
 *  date:   22 February 2023
 * */

#include "benchmark/statbench/analytical_dist.h"

namespace qrc {
namespace benchmark {
namespace statbench {

AnalyticalDistribution::AnalyticalDistribution(const AnalyticalModelData& model_data)
    :ErrorDistribution(model_data.mean, model_data.variance, false),
    model(model_data.model),
    n_events(model_data.n_events),
    error_rate(model_data.error_rate),
    log_error_rate(log(model_data.error_rate)),
    log_error_rate_complement(log(1.0-model_data.error_rate)),
    log_mean(log(model_data.mean))
{}

fp_t
AnalyticalDistribution::pmf(fp_t hw, bool logscale) {
    fp_t hwdiv2 = hw * 0.5;

    fp_t log_prob = 1;
    if (model == AnalyticalModel::binomial) {
        fp_t log_nCk = lgamma(n_events+1) - lgamma(hwdiv2+1) - lgamma(n_events-hwdiv2+1);
        log_prob = log_nCk + hwdiv2 * log_error_rate + (n_events - hwdiv2) * log_error_rate_complement;
    } else if (model == AnalyticalModel::poisson) {
        log_prob = (hwdiv2 - 1) * log_mean - lgamma(hwdiv2+1);
    }

    if (logscale) {
        return log_prob;
    } else {
        return pow(M_E, log_prob);
    }
}

fp_t
AnalyticalDistribution::cdf(fp_t hw, bool logscale) {
    return logscale ? 1.0 : 0.0;
}

AnalyticalDistribution
build_analytical_distribution(const stim::Circuit& circuit, AnalyticalModel model) {
    fp_t sum_error_rate = 0.0;
    fp_t n_error_sources = 0.0;

    DecodingGraph graph = to_decoding_graph(circuit);

    for (auto v : graph.vertices()) {
        for (auto w : graph.adjacency_list(v)) {
            auto e = graph.get_edge(v, w);
            sum_error_rate += e->error_probability;
            n_error_sources += 1.0;
        }
    }
    sum_error_rate *= 0.5;
    n_error_sources *= 0.5;
    fp_t mean_error_rate = sum_error_rate / n_error_sources;

    fp_t mean = sum_error_rate;
    fp_t variance = mean * (1.0 - mean_error_rate);

    AnalyticalModelData d = {
        mean,
        variance,
        n_error_sources,
        mean_error_rate,
        model
    };

    AnalyticalDistribution dist(d);
    return dist;
}

AnalyticalDistribution
bootstrap_from_normal_approximation(const stim::Circuit& circuit, AnalyticalModel model,
    const NumericalDistribution& normal_dist) 
{
    fp_t n_error_sources = 0.0;

    DecodingGraph graph = to_decoding_graph(circuit);

    for (auto v : graph.vertices()) {
        for (auto w : graph.adjacency_list(v)) {
            n_error_sources += 1.0;
        }
    }
    n_error_sources *= 0.5;  // We double counted.
    fp_t mean_error_rate = normal_dist.get_mean() / n_error_sources;
    fp_t mean = normal_dist.get_mean();
    fp_t variance = normal_dist.get_variance();

    AnalyticalModelData d = {
        mean,
        variance,
        n_error_sources,
        mean_error_rate,
        model
    };

    AnalyticalDistribution dist(d);
    return dist;
}

}   // statbench
}   // benchmark
}   // qrc
