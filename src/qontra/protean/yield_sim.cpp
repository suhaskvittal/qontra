/*
 *  author: Suhas Vittal
 *  date:   23 March 2024
 * */

#include "qontra/protean/yield_sim.h"

#include <deque>
#include <set>

#include <math.h>

namespace qontra {
namespace protean {

using namespace net;
using namespace graph;

YieldSimulator::YieldSimulator(PhysicalNetwork* net)
    :freq_map(),
    network(net),
    vcenter(nullptr)
{
    // Compute vertex in the "center" of the processor (no exact center because
    // the processor is 3D). We can use first-degree (# of neighbors) + second-degree
    // (# of neighbors of neighbors) as an alternative.
    size_t max_score = 0;
    for (sptr<phys_vertex_t> v : network->get_vertices()) {
        size_t s = compute_center_score(v);
        if (s > max_score) {
            vcenter = v;
            max_score = s;
        }
    }
}

fp_t
YieldSimulator::est_mean_collisions(
        fp_t prec, uint64_t trials, uint64_t seed, sptr<phys_vertex_t> v)
{
    std::vector<sptr<phys_vertex_t>> vertices;
    if (v == nullptr) {
        vertices = network->get_vertices();
    } else {
        std::set<sptr<phys_vertex_t>> _vertices{v};
        for (sptr<phys_vertex_t> w : network->get_neighbors(v)) {
            _vertices.insert(w);
            for (sptr<phys_vertex_t> u : network->get_neighbors(w)) {
                _vertices.insert(u);
            }
        }
        vertices = std::vector<sptr<phys_vertex_t>>(_vertices.begin(), _vertices.end());
    }

    fp_t mcoll = 0.0;

    std::mt19937_64 rng(seed);
    std::normal_distribution noise{0.0, prec};
    for (uint64_t i = 0; i < trials; i++) {
        std::map<sptr<phys_vertex_t>, fp_t> fmap;
        // Update the frequency map with noise.
        for (auto& [v, f] : freq_map) {
            fmap[v] = freq_map[v] + noise(rng);
        }
        mcoll += static_cast<fp_t>(count_collisions(fmap, vertices));
    }
    mcoll /= static_cast<fp_t>(trials);
    return mcoll;
}

void
YieldSimulator::assign(fp_t prec, const std::vector<fp_t>& freq_list) {
    std::cout << "[ yield_sim ]\n";
    std::deque<sptr<phys_vertex_t>> bfs{vcenter};
    std::set<sptr<phys_vertex_t>> visited;
    while (bfs.size()) {
        sptr<phys_vertex_t> v = bfs.front();
        bfs.pop_front();
        if (visited.count(v)) {
            continue;
        }
        if (freq_map.empty()) {
            // Set frequency to be middle of freq_list.
            freq_map[v] = freq_list.at(freq_list.size() / 2);
        } else {
            fp_t min_coll = std::numeric_limits<fp_t>::max();
            fp_t best_f = -1;
            for (fp_t f : freq_list) {
                // Simulate the yield of the processor locally.
                freq_map[v] = f;
                fp_t c = est_mean_collisions(prec, 1000, 0, v);
                if (c < min_coll) {
                    min_coll = c;
                    best_f = f;
                }
            }
            if (best_f < 0) {
                std::cerr << "no good frequency found for " << print_v(v) << std::endl;
                exit(1);
            }
            freq_map[v] = best_f;
        }
        std::cout << "\tassigned F=" << freq_map[v] << " ---> " << print_v(v) << std::endl;
        for (sptr<phys_vertex_t> w : network->get_neighbors(v)) {
            bfs.push_back(w);
        }
        visited.insert(v);
    }
}

size_t
YieldSimulator::count_violations(
        sptr<phys_vertex_t> v,
        const std::map<sptr<phys_vertex_t>, fp_t>& fmap) 
{
    size_t violations = 0;
    const auto adj = network->get_neighbors(v);
    for (sptr<phys_vertex_t> w : adj) {
        if (!fmap.count(v) || !fmap.count(w)) continue;
        // Conditions:
        const bool c1 = abs(fmap.at(v) - fmap.at(w)) <= 17e6,
                   c2 = abs(fmap.at(v) - (fmap.at(w)-0.5*ANH)) <= 4e6,
                   c3 = abs(fmap.at(v) - (fmap.at(w)-ANH)) <= 25e6,
                   c4 = fmap.at(v) > fmap.at(w)-ANH;
        violations += static_cast<size_t>(c1) 
                        + static_cast<size_t>(c2)
                        + static_cast<size_t>(c3)
                        + static_cast<size_t>(c4);
    }
    for (size_t i = 0; i < adj.size(); i++) {
        sptr<phys_vertex_t> w = adj.at(i);
        if (!fmap.count(w)) continue;
        for (size_t j = i+1; j < adj.size(); j++) {
            sptr<phys_vertex_t> u = adj.at(j);
            if (!fmap.count(u)) continue;
            const bool c4 = abs(fmap.at(w) - fmap.at(u)) <= 17e6,
                       c5 = abs(fmap.at(w) - (fmap.at(u)-ANH)) <= 25e6;
            const bool c6 = fmap.count(v)
                            && abs((2*fmap.at(v)+ANH) - (fmap.at(w)+fmap.at(u))) <= 17e6;
            violations += 
                static_cast<size_t>(c4) + static_cast<size_t>(c5) + static_cast<size_t>(c6);
        }
    }
    return (violations > 0);
}

size_t
YieldSimulator::count_collisions(
        const std::map<sptr<phys_vertex_t>, fp_t>& fmap,
        const std::vector<sptr<phys_vertex_t>>& vertices)
{
    size_t coll = 0;
    for (sptr<phys_vertex_t> v : vertices) {
        coll += count_violations(v, fmap);
    }
    return coll;
}

size_t
YieldSimulator::compute_center_score(sptr<phys_vertex_t> v) {
    size_t score = 0;
    for (sptr<phys_vertex_t> w : network->get_neighbors(v)) {
        score++;
        for (sptr<phys_vertex_t> u : network->get_neighbors(w)) {
            score += (u != v && !network->contains(u, v));
        }
    }
    return score;
}

}   // protean
}   // qontra
