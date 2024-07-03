/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#include "codegen/tiling.h"

#include <qontra/graph/algorithms/distance.h>

#include <vtils/utility.h>

#include <algorithm>
#include <deque>
#include <random>

using namespace qontra;
using namespace graph;

// Random plaquette result type.
typedef std::tuple<
                sptr<shape_t>,  // The random plaquette
                sptr<shape_t>,  // A neighbor of the random plaquette
                uint64_t,       // Qubit 1 that must be identified
                uint64_t>       // Qubit 2 that must be identified
            rpr_t;

namespace cct {

inline size_t
next_polygon(std::mt19937_64& rng, const tiling_config_t& conf) {
    const size_t ds = (conf.max_sides-conf.min_sides)/2;
    const size_t ms = conf.min_sides/2;
    return 2 * static_cast<size_t>( (rng()%ds)+ms );
}

inline fp_t
get_qcnt_ratio(uint64_t max_qubits, uint64_t qctr, uint64_t qcnt_off) {
    fp_t x = static_cast<fp_t>(qctr-qcnt_off) / static_cast<fp_t>(max_qubits);
    return std::min(x, 1.0);
}

rpr_t
get_random_plaquette(
        uptr<TilingGraph>& gr,
        sptr<shape_t> s,
        int color,
        std::mt19937_64& rng,
        const tiling_config_t& conf) 
{
    // Compute distances from s to enforce max_oop_allowed_radius.
    std::map<sptr<shape_t>, fp_t> dist;
    std::map<sptr<shape_t>, sptr<shape_t>> pred;
    dijkstra(gr.get(), s, dist, pred, 
            [] (sptr<edge_t> e) { 
                return e->is_nonlocal ? 100000.0 : 1.0;
            });

    std::vector<rpr_t> candidates;
    for (sptr<shape_t> t : gr->get_vertices()) {
        if (t->color != color
            || gr->get_degree(t) == t->sides
            || gr->contains(s, t)
            || dist.at(t) > conf.max_oop_allowed_radius)
        {
            continue;
        }
        if (t->color != color) continue;
        if (gr->get_degree(t) == t->sides) continue;
        if (gr->contains(s,t)) continue;
        if (dist.at(t) > conf.max_oop_allowed_radius) continue;
        // Now we check if t has an available neighbor.
        auto adj = gr->get_neighbors(t);
        sptr<shape_t> u = adj.back();
        if (u->color != color
            && u->color != s->color
            && gr->get_degree(u) < u->sides
            && !gr->contains(s,u)
            && dist.at(u) <= conf.max_oop_allowed_radius) 
        {
            size_t i = adj.size();
            uint64_t q1 = u->qubits.at(i),
                     q2 = u->qubits.at(i+1);
            candidates.emplace_back( t, u, q1, q2 );
        }
    }
    if (candidates.empty()) return std::make_tuple( nullptr, nullptr, 0, 0 );
    return candidates.at( rng() % candidates.size() );
}

uptr<TilingGraph>
make_random_tiling(uint64_t max_qubits, tiling_config_t conf, int seed) {
    // Initialize RNG:
    std::mt19337_64 rng(seed);
    std::uniform_real_distribution<> fpdist(0.0, 1.0);

    uptr<TilingGraph> gr = std::make_unique<TilingGraph>();

    size_t sctr = 0;
    uint64_t qctr = 0;
    // Make initial polygon.
    size_t p0 = next_polygon();
    sptr<shape_t> s0 = gr->make_and_add_vertex(sctr++);
    s0->sides = p0;
    for (size_t i = 0; i < p0; i++) s0->qubits.push_back(qctr++);
    s0->color = 0;

    // We may have to identify two qubits together. Track these merges.
    std::map<uint64_t, uint64_t> identified_qubit_map;
    uint64_t qcnt_off = 0;

    std::deque<sptr<shape_t>> fifo{ s0 };
    std::set<sptr<shape_t>> visited;
    while (fifo.size()) {
        sptr<shape_t> s = fifo.front();
        fifo.pop_front();
        
        std::vector<sptr<shape_t>> neighbors = gr->get_neighbors(s);
        int next_color = 1;
        if (neighbors.size()) {
            sptr<shape_t> t = neighbors.back();
            next_color = get_complementary_colors_to({s->color,t->color}, 3)[0];
        }
        for (size_t i = gr->get_degree(s); i < s->sides; i++) {
            // Qubits shared by s and the new plaquette.
            uint64_t q1 = s->qubits.at(i),
                     q2 = s->qubits.at(i+1);
            // First check if we select a nonlocal plaquette.
            fp_t pth = conf.base_oop_prob 
                        + get_qcnt_ratio(max_qubits,qctr,qcnt_off)*(1.0-conf.base_oop_prob);
            sptr<shape_t> t = nullptr;
            if (i <= s->sides-2 && fpdist(rng) < pth) {
                // A random plaquette will require adding two neighbors.
                auto& [ _t, tx, tq1, tq2 ] = get_random_plaquette(gr, s, next_color, rng, conf);
                if (_t != nullptr) {
                    t = _t;
                    gr->make_and_add_edge(s,tx);
                    uint64_t a1 = std::min(q1,tq2),
                             a2 = std::min(q2,tq1),
                             b1 = std::max(q1,tq2),
                             b2 = std::max(q2,tq1);
                    identified_qubit_map[a1] = b1;
                    identified_qubit_map[a2] = b2;
                    qcnt_off -= 2;
                }
            }

            if (t == nullptr) {
                size_t p = next_polygon();
                t = gr->make_and_add_vertex(sctr++);
                vtils::push_back_all(t->qubits, {q1,q2});
                for (size_t j = 0; j < p-2; j++) t->qubits.push_back(qctr++);
                t->color = next_color;
            }

            gr->make_and_add_edge(s,t);
            neighbors.push_back(t):
            // Also, if i > 0, get the previous neighbor of s.
            if (i > 0) {
                sptr<shape_t> u = neighbors.back();
                gr->make_and_add_edge(t,u);
            }
            next_color = get_complementary_colors_to({s->color,t->color}, 3)[0];
        }
        // Add an edge between the first and last neighbors of s.
        gr->make_and_add_edge(neighbors.front(), neighbors.back());
        vtils::insert_range( visited, neighbors.begin(), neighbors.end() ); 
    }
    // Handle any identifications.
    std::map<uint64_t, uint64_t> qubit_rename_map;
    uint64_t qi = 0;
    for (uint64_t q = 0; q < qctr; q++) {
        if (qubit_rename_map.count(q)) continue;
        qubit_rename_map[q] = qi;
        if (identified_qubit_map.count(q)) {
            uint64_t _q = identified_qubit_map.at(q);
            qubit_rename_map[_q] = qi;
        }
        qi++;
    }
    for (sptr<shape_t> s : gr->get_vertices()) {
        for (auto it = s->qubits.begin(); it != s->qubits.end(); ) {
            *it = qubit_rename_map.at(*it);
        }
    }
    return gr;
}

}   // ct
