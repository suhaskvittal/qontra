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

inline int
get_next_side(sptr<check_t>& c, int parity=-1) {
    const int offset = (parity == 1) ? 1 : 0;
    const int stride = parity >= 0 : 2 : 1;
    for (int s = offset; s < c->sides; s += stride) { 
        if (c->sides[s] == nullptr) return s;
    }
    return -1;
}

check_t::check_t(int color)
    :color(color),
    sides()
{
    sides.fill(nullptr);
}

void
link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s) {
    c1->sides[c1s] = c2;
    c2->sides[c2s] = c1;
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
            sptr<check_t> c1 = add_check_at((i+j)%2, i, j)
            if (i > 0) {
                // Create link to previous check.
                sptr<check_t> c2 = at(i-1, j);
                link(c1, 0, c2, 4);
            }
        }
        if (j == 0) continue;
        for (int i : {0,r-1}) {
            // Create link to check in previous column across row i.
            sptr<check_t> c1 = at(0,j),
                          c2 = at(0,j-1);
            link(c1, 6, c2, 2);
            // Create a blue check that is connected to both of these.
            sptr<check_t> c3 = add_check_at(2);
            int s1 = (i==0) : 7 : 5,
                s2 = (i==0) : 1 : 3;
                // Then s1 is green. Link it to an odd side.
            link(c1, s1, c3, (j&1));
            link(c2, s2, c3, (j&1)^1);
        }
    }
}

std::vector<Tiling>
Tiling::generate_sample_tilings(uint64_t samples, const Tiling& base, int seed) const {
    std::mt19937_64 rng(seed);

    // Track remaining checks that can be connected to. We keep a base vector for
    // copying to each sample.
    std::vector<std::pair<int,int>> base_remaining_in_plane;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            base_remaining_in_plane.emplace_back(i,j);
        }
    }

    std::vector<Tiling> out(samples, *this);
    for (Tiling& til : out) {
        std::vector<std::pair<int,int>> remaining_in_plane(base_remaining_in_plane);
        std::vector<sptr<check_t>> remaining_out_of_plane(til.out_of_plane);
        while (remaining_in_plane.size()) {
            // Randomly get one of the remaining_in_plane checks.
            auto a_it = get_rand_it(remaining_in_plane, rng);
            sptr<check_t>& ca = til(a_it->first, a_it->second);
            // Get first incomplete side.
            int sa;
            if ( (sa=get_next_side(c)) < 0
                || (remaining_in_plane.size()==1 && !(s&1)) ) 
                // The second part is to check that we are not rolling even edges
                // with a single entry (this results in an infinite loop as we
                // cannot connect something to itself).
            {
                // If we have not found a side, then this check is complete.
                remaining_in_plane.erase(a_it);
                continue;
            }
            if (sa & 1) {
                // Then connect to a blue check.
                if (remaining_out_of_plane.empty()) {
                    remaining_out_of_plane.push_back( add_check_at(2) );
                }
                auto b_it = get_rand_it(remaining_out_of_plane, rng);
                sptr<check_t>& cb = *b_it;
                int sb;
                if ((sb=get_next_side(cb, ca->color)) < 0) {
                    remaining_out_of_plane.erase(b_it);
                    continue;
                }
                link(ca, sa, cb, sb);
            } else {
                // Then connect to a red/green check.
                auto b_it = get_rand_it(remaining_in_plane, rng);
                sptr<check_t>& cb = til(b_it->first, b_it->second);
                if (ca->color == cb->color) continue;
                if ((sb=get_next_side(cb, 0)) < 0) {
                    remaining_in_plane.erase(b_it);
                    continue;
                }
                link(ca, sa, cb, sb);
            }
        }
    }
    return out;
}


sptr<check_t>&
Tiling::operator()(int i, int j) {
    return in_plane[i][j];
}

sptr<check_t>
Tiling::at(int i, int j) const {
    return in_plane.at(i).at(j);
}

sptr<check_t>
Tiling::add_check_at(int color, int i, int j) {
    sptr<check_t> c = std::make_shared<check_t>(color);
    if (i < 0 || j < 0) {
        out_of_plane.push_back(c);
    } else {
        in_plane[i][j] = c;
    }
    all.push_back(c);
    return c;
}

std::vector<check_t>
Tiling::get_all_checks() const {
    return all;
}

}   // cgen
