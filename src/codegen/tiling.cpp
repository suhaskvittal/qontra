/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

#include "codegen/tiling.h"

#include <random>

namespace cgen {

template <class CONTAINER, ITER> inline ITER
get_rand_it(const CONTAINER& x, std::mt19937_64& rng) {
    return x.begin() + (rng() % x.size());
}


check_t::check_t(int color, int s, int i, int j)
    :color(color),
    sides(s, nullptr),
    i(i),
    j(j)
{}

int
get_side(sptr<check_t> a, sptr<check_t> b, int offset, bool next_must_be_empty=false) {
    for (int i = offset; i < a->size(); i += 2) {
        if (a->get(i) == nullptr
            && a->get(i-2) != b
            && a->get(i+2) != b
            && (!next_must_be_empty || a->get(i+1) == nullptr)) 
        {
            return i;
        }
    }
    return -1;
}

Tiling::Tiling(int r, int c)
    :r(r),
    c(c),
    in_plane( r, std::vector<sptr<check_t>>(c, nullptr) ),
    out_of_plane(),
    all()
{
    for (int j = 0; j < c; j++) {
        for (int i = 0; i < r; i++) {
            sptr<check_t> c1 = add_check_at((i+j)%2, 8, i, j)
            if (i > 0) {
                // Create link to previous check.
                sptr<check_t> c2 = at(i-1, j);
                link(c1, 0, c2, 4);
            }
        }
    }
}

std::vector<Tiling>
Tiling::generate_sample_tilings(uint64_t samples, const Tiling& base, int seed) const {
    std::mt19937_64 rng(seed);

    std::vector<Tiling> out(samples, *this);
    for (Tiling& til : out) {
        std::vector<sptr<check_t>> remaining(base.in_place);
        while (remaining.size()) {
            sptr<check_t> c_oop = add_check_at(2, blue_sides);
            sptr<check_t> c_init;
            int s_init;
            do {
                auto it = get_rand_it(remaining, rng);
                c_init = at(it->first, it->second);
                s_init = get_side(c_init, c_oop, 1);
            } while (s_init < 0);
            int s_oop = c_init->color;
            link(c_oop, s_oop, c_init, s_init);

            sptr<check_t> c_curr = c_init;
            int s_curr = s_init;
            for (int del = 1; del < c_oop->sides.size(); del++) {
                if (c_curr->get(s_curr+1) == nullptr) {
                    // Find a neighbor for c_curr.
find_neighbor:
                    auto it = get_rand_it(remaining, rng);
                    if (*it == c_curr
                            || c_curr->get(s_curr-1) == *it
                            || c_curr->get(s_curr+3) == *it)
                    {
                        goto find_neighbor;
                    }
                    int s = get_side(*it, c_curr, 0, true);
                    if (s < 0) goto find_neighbor;
                    // Otherwise, this neighbor is fine.
                    link(c_curr, s_curr+1, *it, s);
                    link(c_oop, s_oop+del, *it, s+1);
                    // Update curr.
                    c_curr = *it;
                    s_curr = s+1;
                } else {
                    sptr<check_t> next = c_curr->get(s_curr+1);
                    // This will get the corresponding side on which c_curr lies.
                    int s = next->side_map[c_curr][s_curr+1];
                    link(c_oop, s_oop+del, next, s+1);
                    c_curr = next;
                    s_curr = s+1;
                }
            }
            // Now, we must connect c_curr to c_init.
        }
    }
    return out;
}

vtils::Mat2
Tiling::to_matrix() const {
    // Identify data qubits.
    int n = 0;
    std::map<sptr<check_t>, std::vector<int>> qubit_map;
    for (const sptr<check_t>& c : all) {
        qubit_map[c] = std::vector<int>(8, -1);
    }
    for (const sptr<check_t>& c : all) {
        
    }
}

}   // cgen
