/*
 *  author: Suhas Vittal
 *  date:   2 July 2024
 * */

#include "codegen/tiling.h"

#include <qontra/graph/algorithms/distance.h>
#include <qontra/graph/decoding_graph.h>

#include <vtils/utility.h>

#include <algorithm>
#include <deque>
#include <random>

using namespace qontra;
using namespace graph;
using namespace vtils;

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

inline void
identify(std::map<uint64_t, uint64_t>& m, uint64_t q1, uint64_t q2) {
    if (q1 > q2) std::swap(q1,q2);
    m[q1] = q2;
}

inline sptr<shape_t>
make_and_add_shape(uptr<TilingGraph>& g, uint64_t& sctr, size_t p, int c) {
    sptr<shape_t> s = g->make_and_add_vertex(sctr++);
    s->sides = p;
    s->qubits = std::vector<uint64_t>(p);
    s->neighbors = std::vector<sptr<shape_t>>(p, nullptr);
    s->fptr = 0;
    s->bptr = p-1;
    s->color = c;
    return s;
}

inline sptr<edge_t>
make_and_add_edge(uptr<TilingGraph>& g, sptr<shape_t> s, sptr<shape_t> t, bool sf, bool tf) {
    if (sf) s->set_neighbor(s->fptr++, t);
    else    s->set_neighbor(s->bptr--, t);
    if (tf) t->set_neighbor(t->fptr++, s);
    else    t->set_neighbor(t->bptr--, s);
    return g->make_and_add_edge(s,t);
}

inline void
identify_common_in_triangle(
        std::map<uint64_t, uint64_t>& identified_qubit_map,
        sptr<shape_t> s, sptr<shape_t> t, sptr<shape_t> u,
        uint64_t q1, uint64_t q2, uint64_t& qcnt_off) 
{
    auto t_it = std::find(t->qubits.begin(), t->qubits.end(), q1);
    auto u_it = std::find(u->qubits.begin(), u->qubits.end(), q1);
    if (*(t_it+1) == q2) {
        t_it = t_it == t->qubits.begin() ? t->qubits.end()-1 : t_it-1;
        u_it = u_it == u->qubits.end()-1 ? u->qubits.begin() : u_it+1;
    } else {
        t_it = t_it == t->qubits.end()-1 ? t->qubits.begin() : t_it+1;
        u_it = u_it == u->qubits.begin() ? u->qubits.end()-1 : u_it-1;
    }
    identify(identified_qubit_map, *t_it, *u_it);
    std::cout << "(" << q1 << "," << q2 << ") identified " 
        << *t_it << " <---> " << *u_it << " for "
        << print_v(s) << "," << print_v(t) << print_v(u) << std::endl;
    qcnt_off -= 2;
}

inline void
identify_common_in_nonlocal_triangle(
        std::map<uint64_t, uint64_t>& identified_qubit_map,
        sptr<shape_t> s, sptr<shape_t> t, sptr<shape_t> u,
        uint64_t q1, uint64_t q2, uint64_t tq1, uint64_t tq2, 
        uint64_t& qcnt_off) 
{
    auto t_it = std::find(t->qubits.begin(), t->qubits.end(), tq1);
    auto u_it = std::find(u->qubits.begin(), u->qubits.end(), q1);
    if (*(t_it+1) == tq2) {
        t_it = t_it == t->qubits.begin() ? t->qubits.end()-1 : t_it-1;
        u_it = u_it == u->qubits.end()-1 ? u->qubits.begin() : u_it+1;
    } else {
        t_it = t_it == t->qubits.end()-1 ? t->qubits.begin() : t_it+1;
        u_it = u_it == u->qubits.begin() ? u->qubits.end()-1 : u_it-1;
    }
    identify(identified_qubit_map, *t_it, *u_it);
    std::cout << "NL(" << q1 << "," << q2 << ") identified " 
        << *t_it << " <---> " << *u_it << " for "
        << print_v(s) << "," << print_v(t) << print_v(u) << std::endl;
    qcnt_off -= 2;
}

void
rename(
    const std::map<uint64_t, uint64_t>& identified_qubit_map,
    std::vector<sptr<shape_t>> vertices,
    size_t qctr) 
{
    // Handle any identifications.
    std::map<uint64_t, uint64_t> qubit_rename_map;
    uint64_t qi = 0;
    for (uint64_t q = 0; q < qctr; q++) {
        if (qubit_rename_map.count(q)) continue;
        qubit_rename_map[q] = qi;
        uint64_t _q = q;
        while (identified_qubit_map.count(_q)) {
            _q = identified_qubit_map.at(_q);
            qubit_rename_map[_q] = qi;
        }
        qi++;
    }
    for (sptr<shape_t> s : vertices) {
        for (auto it = s->qubits.begin(); it != s->qubits.end(); it++) {
            *it = qubit_rename_map.at(*it);
        }
    }
}

sptr<shape_t>
select_random_plaquette(
        uptr<TilingGraph>& g, sptr<shape_t> s, int color,
        std::mt19937_64& rng, const tiling_config_t& conf)
{
    std::vector<sptr<shape_t>> candidates;
    sptr<shape_t> t = s->get_neighbor(s->fptr-1);
    for (sptr<shape_t> x : g->get_vertices()) {
        // TODO: test locality.
        if (x->color != color
            || g->get_degree(x) > x->sides-3
            || g->contains(x,s))
        {
            continue;
        }
        sptr<shape_t> y = s->get_neighbor(x->bptr+1);

        candidates.push_back(x);
    }
    if (candidates.empty()) return nullptr;
    return candidates[rng() % candidates.size()];
}

uptr<TilingGraph>
make_random_tiling(uint64_t max_qubits, tiling_config_t conf, int seed) {
    // Initialize RNG:
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> fpdist(0.0, 1.0);

    uptr<TilingGraph> gr = std::make_unique<TilingGraph>();

    uint64_t sctr = 0,
             qctr = 0;
    // Make initial polygon.
    sptr<shape_t> s0 = make_and_add_shape(gr, sctr, next_polygon(rng, conf), 0);
    for (size_t i = 0; i < s0->sides; i++) {
        s0->set_qubit(i, qctr++);
    }

    // We may have to identify two qubits together. Track these merges.
    std::map<uint64_t, uint64_t> identified_qubit_map;
    uint64_t qcnt_off = 0;

    std::deque<sptr<shape_t>> fifo{ s0 };
    std::set<sptr<shape_t>> visited;
    while (fifo.size()) {
        sptr<shape_t> s = fifo.front();
        fifo.pop_front();
        if (visited.count(s)) continue;

        std::vector<int> color_array(s->sides);
        if (gr->get_degree(s) == 0) {
            std::vector<int> cc = get_complementary_colors_to({s->color}, 3);
            for (size_t i = 0; i < s->sides; i++) {
                color_array[i] = cc.at(i&1);
            }
            // Also make an initial neighboring plaquette.
            sptr<shape_t> t0 = 
                make_and_add_shape(gr, sctr, next_polygon(rng,conf), color_array[0]);
            for (size_t i = 0; i < t0->sides; i++) {
                if (i <= 1) t0->set_qubit(i, s->get_qubit(i));
                else        t0->set_qubit(i, qctr++);
            }
            make_and_add_edge(gr, s, t0, true, true);
        } else {
            int c = s->get_neighbor(s->fptr-1)->color;
            color_array[s->fptr-1] = c;
            for (size_t i = 0; i < s->sides; i++) {
                int x = i&1 ? c : get_complementary_colors_to({s->color,c},3)[0];
                color_array[(s->fptr+i-1) % s->sides] = x;
            } 
        }
        // Now that we know the colors in the plaquettes, it's time to make the
        // neighboring plaquettes.
        bool last_was_nonlocal = false;
        size_t init_fptr = s->fptr;
        for (size_t i = init_fptr; i < init_fptr+s->sides; i++) {
            int c = color_array.at(i % s->sides);
            sptr<shape_t> t = s->get_neighbor(i);
            sptr<shape_t> u = s->get_neighbor(i-1);

            uint64_t q1 = s->get_qubit(i),
                     q2 = s->get_qubit(i+1);
    
            rename(identified_qubit_map, gr->get_vertices(), qctr);
            identified_qubit_map.clear();

            if (t == nullptr) {
                // Select whether or not we want to use an existing plaquette or
                // make a new one.
                fp_t pth = conf.base_oop_prob 
                            + (1.0-conf.base_oop_prob)
                                *get_qcnt_ratio(max_qubits, qctr, qcnt_off);
                if (!last_was_nonlocal && fpdist(rng) < pth
                        && (t=select_random_plaquette(gr,s,c,rng,conf)) != nullptr) 
                {
                    // Identify qubits between s and t.
                    uint64_t tq1 = t->get_qubit(t->bptr-1),
                             tq2 = t->get_qubit(t->bptr-2);
                    identify(identified_qubit_map, q1, tq1);
                    identify(identified_qubit_map, q2, tq2);

                    // u gets connected to v.
                    make_and_add_edge(gr, s, t, true, false);
                    make_and_add_edge(gr, t, u, false, !last_was_nonlocal);

                    identify_common_in_nonlocal_triangle(
                            identified_qubit_map,s,t,u,q1,q2,tq1,tq2,qcnt_off);
                    last_was_nonlocal = true;
                } else {
                    t = make_and_add_shape(gr, sctr, next_polygon(rng,conf), c);
                    t->set_qubit(0,q1);
                    t->set_qubit(1,q2);
                    for (size_t j = 2; j < t->sides; j++) {
                        t->set_qubit(j,qctr++);
                    }
                    make_and_add_edge(gr, s, t, true, true);
                    make_and_add_edge(gr, t, u, false, !last_was_nonlocal);

                    identify_common_in_triangle(identified_qubit_map,s,t,u,q1,q2,qcnt_off);
                    last_was_nonlocal = false;
                }
            } else if (!gr->contains(t,u)) {
                make_and_add_edge(gr, t, u, false, !last_was_nonlocal);
                identify_common_in_triangle(identified_qubit_map,s,t,u,q1,q2,qcnt_off);
                last_was_nonlocal = false;
            }
            fifo.push_back(t);
        }
        visited.insert(s);
    }
    return gr;
}

}   // ct
