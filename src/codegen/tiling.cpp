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

typedef std::tuple<sptr<shape_t>, sptr<shape_t>, sptr<shape_t>> face_t;

inline size_t
next_polygon(std::mt19937_64& rng, const tiling_config_t& conf) {
    if (conf.max_sides == conf.min_sides) return conf.max_sides;
    const size_t ds = (conf.max_sides-conf.min_sides)/2;
    const size_t ms = conf.min_sides/2;
    return 2 * static_cast<size_t>( (rng()%ds)+ms );
}

inline fp_t
get_qcnt_ratio(uint64_t max_qubits, uint64_t q) {
    fp_t x = static_cast<fp_t>(q) / static_cast<fp_t>(max_qubits);
    return std::min(x, 1.0);
}

inline sptr<shape_t>
make_and_add_shape(uptr<TilingGraph>& g, uint64_t& sctr, size_t p, int c) {
    sptr<shape_t> s = g->make_and_add_vertex(sctr++);
    s->sides = p;
    s->qubits = std::vector<uint64_t>(p);
    s->neighbors = std::vector<sptr<shape_t>>(p, nullptr);
    s->cnt = 0;
    s->fptr = 0;
    s->bptr = p-1;
    s->color = c;
    return s;
}

inline sptr<edge_t>
make_and_add_edge(uptr<TilingGraph>& g, sptr<shape_t> s, sptr<shape_t> t, bool sf, bool tf) {
    if (s->cnt == s->sides) {
        std::cerr << "s(" << print_v(s) << ") degree overflow error." << std::endl;
        exit(1);
    }
    if (t->cnt == t->sides) {
        std::cerr << "t(" << print_v(t) << ") degree overflow error." << std::endl;
        exit(1);
    }
    if (s->get_neighbor(s->fptr) != nullptr) {
        std::cerr << "s(" << print_v(s) << ") neighbor already exists." << std::endl;
        exit(1);
    }
    if (t->get_neighbor(t->fptr) != nullptr) {
        std::cerr << "t(" << print_v(t) << ") neighbor already exists." << std::endl;
        exit(1);
    }
    if (sf) s->set_neighbor(s->fptr++, t);
    else    s->set_neighbor(s->bptr--, t);
    if (tf) t->set_neighbor(t->fptr++, s);
    else    t->set_neighbor(t->bptr--, s);
    return g->make_and_add_edge(s,t);
}

inline face_t
make_face(sptr<shape_t> x, sptr<shape_t> y, sptr<shape_t> z) {
    if (x < y && y < z)         return std::make_tuple(x,y,z);
    else if (x < z && z < y)    return std::make_tuple(x,z,y);
    else if (y < x && x < z)    return std::make_tuple(y,x,z);
    else if (y < z && z < x)    return std::make_tuple(y,z,x);
    else if (z < x && x < y)    return std::make_tuple(z,x,y);
    else                        return std::make_tuple(z,y,x);
}

sptr<shape_t>
select_random_plaquette(
        uptr<TilingGraph>& g, sptr<shape_t> s, int color,
        std::mt19937_64& rng, const tiling_config_t& conf)
{
    std::vector<sptr<shape_t>> candidates;
    sptr<shape_t> t = s->get_neighbor(s->fptr-1);

    if (s->cnt > s->sides - 3) return nullptr;
    if (t->cnt > t->sides - 2) return nullptr;

    for (sptr<shape_t> x : g->get_vertices()) {
        // TODO: test locality.
        if (x->color != color
            || x->cnt > x->sides-3
            || (g->contains(x,s) && g->contains(x,t)))
        {
            continue;
        }
        sptr<shape_t> y = s->get_neighbor(x->bptr+1);

        candidates.push_back(x);
    }
    if (candidates.empty()) return nullptr;
    return candidates[rng() % candidates.size()];
}

void
assign_qubits(uptr<TilingGraph>& g) {
    std::map<face_t, uint64_t> face_qubit_map;
    uint64_t qctr = 0;
    auto vertices = g->get_vertices();
    for (sptr<shape_t> s : vertices) {
        for (size_t i = 0; i < s->sides; i++) {
            sptr<shape_t> t = s->get_neighbor(i);
            sptr<shape_t> u = s->get_neighbor(i+1);
            face_t f = make_face(s,t,u);
            if (t == nullptr || u == nullptr) {
                std::cerr << "found null element in face" << std::endl;
                exit(1);
            }

            if (!face_qubit_map.count(f)) {
                std::cout << "face " << print_v(std::get<0>(f)) << "|" 
                    << print_v(std::get<1>(f)) << "|" 
                    << print_v(std::get<2>(f)) << " --> " << qctr << std::endl;
                face_qubit_map[f] = qctr++;
            }
            s->qubits[i] = face_qubit_map[f];
        }
    }
}

uptr<TilingGraph>
make_random_tiling(uint64_t max_qubits, tiling_config_t conf, int seed) {
    // Initialize RNG:
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> fpdist(0.0, 1.0);

    uptr<TilingGraph> gr = std::make_unique<TilingGraph>();

    uint64_t sctr = 0;
    // Make initial polygon.
    sptr<shape_t> s0 = make_and_add_shape(gr, sctr, next_polygon(rng, conf), 0);

    std::deque<sptr<shape_t>> fifo{ s0 };
    std::set<sptr<shape_t>> visited;
    while (fifo.size()) {
        sptr<shape_t> s = fifo.front();
        fifo.pop_front();
        if (visited.count(s)) continue;

        std::vector<int> color_array(s->sides);
        if (s->cnt == 0) {
            std::vector<int> cc = get_complementary_colors_to({s->color}, 3);
            for (size_t i = 0; i < s->sides; i++) {
                color_array[i] = cc.at(i&1);
            }
            // Also make an initial neighboring plaquette.
            sptr<shape_t> t0 = 
                make_and_add_shape(gr, sctr, next_polygon(rng,conf), color_array[0]);
            make_and_add_edge(gr, s, t0, true, true);
        } else {
            int c = s->get_neighbor(s->fptr-1)->color;
            color_array[s->fptr-1] = c;
            for (size_t i = 1; i < s->sides; i++) {
                int x = i&1 ? get_complementary_colors_to({s->color,c},3)[0] : c;
                color_array[(s->fptr+i-1) % s->sides] = x;
            } 
        }

        // Now that we know the colors in the plaquettes, it's time to make the
        // neighboring plaquettes.
        bool last_was_nonlocal = false;
        size_t init_fptr = s->fptr;

        std::cout << print_v(s) << ":" << std::endl;
        std::cout << "\tPREV:";
        for (size_t i = init_fptr; i < init_fptr+s->sides; i++) {
            std::cout << " " << print_v(s->get_neighbor(i));
        }
        std::cout << std::endl;

        std::cout << "\tNEW:";

        for (size_t i = init_fptr; i < init_fptr+s->sides; i++) {
            int c = color_array.at(i % s->sides);
            sptr<shape_t> t = s->get_neighbor(i);
            sptr<shape_t> u = s->get_neighbor(i-1);

            if (t == nullptr) {
                // Select whether or not we want to use an existing plaquette or
                // make a new one.
                fp_t pth = conf.base_oop_prob 
                            + (1.0-conf.base_oop_prob)
                                *get_qcnt_ratio(max_qubits, gr->m()/2);
                if (fpdist(rng) < pth
                        && (t=select_random_plaquette(gr,s,c,rng,conf)) != nullptr) 
                {
                    std::cout << " " << print_v(t) << "(r)";
                    // For t, we may do something a bit more manual.
                    gr->make_and_add_edge(s,t);
                    gr->make_and_add_edge(t,u);

                    s->set_neighbor(s->fptr++, t);
                    t->set_neighbor(t->bptr+1, s);
                    t->set_neighbor(t->bptr+2, u);
                    if (last_was_nonlocal) {
                        u->set_neighbor(u->bptr--, t);
                    } else {
                        u->set_neighbor(u->fptr++, t);
                    }
                } else {
                    t = make_and_add_shape(gr, sctr, next_polygon(rng,conf), c);
                    std::cout << " " << print_v(t) << "(n)";
                    make_and_add_edge(gr, s, t, true, true);
                    make_and_add_edge(gr, t, u, false, !last_was_nonlocal);
                    last_was_nonlocal = false;
                }
            } else if (!gr->contains(t,u)) {
                std::cout << " " << print_v(t) << "(o)";
                make_and_add_edge(gr, t, u, false, !last_was_nonlocal);
                last_was_nonlocal = false;
            } else {
                std::cout << " " << print_v(t) << "(x)";
            }
            fifo.push_back(t);
        }
        std::cout << std::endl;
        visited.insert(s);
    }
    // Assign qubits to each plaquette.
    assign_qubits(gr);
    return gr;
}

}   // ct
