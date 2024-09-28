/*
 *  author: Suhas Vittal
 *  date:   8 August 2024
 * */
#define DEBUG

#include "codegen/monte_carlo.h"

#include <vtils/utility.h>

#include <algorithm>
#include <deque>
#include <map>
#include <set>

namespace cgen {

inline bool has_free_odd_side(sptr<check_t> c) {
    for (int s = 1; s < c->size(); s += 2) {
        if (c->get(s) == nullptr) return true;
    }
    return false;
}

MonteCarloManager::MonteCarloManager(int seed)
    :rng(seed),
    fp_distr(0.0, 1.0)
{}

Tiling
MonteCarloManager::make_rand_init_tiling() {
    const int& r = config.r;
    const int& c = config.c;
    Tiling t(r, c);

    for (int j = 0; j < c; j++) {
        for (int i = 0; i < r; i++) {
            sptr<check_t> c1 = t.add_check_at((i+j)%2, 8, i, j);
            if (i > 0 && fp_distr(rng) <= config.planar_link_rate) {
                // Create link to previous check.
                sptr<check_t> c2 = t.at(i-1, j);
                link(c1, 0, c2, 4);
            }
        }
        if (j == 0) continue;
        for (int i = 0; i < r; i++) {
            if (fp_distr(rng) <= config.planar_link_rate) {
                // Create link to check in previous column across row i.
                sptr<check_t> c1 = t.at(i,j),
                              c2 = t.at(i,j-1);
                link(c1, 6, c2, 2);
            }
        }
    }
    return t;
}

std::vector<Tiling>
MonteCarloManager::run(uint64_t n) {
    std::vector<Tiling> samples;
    for (uint64_t i = 0; i < n; i++) {
        Tiling t = make_rand_init_tiling();

        bool found_empty = false;
        for (sptr<check_t> c_main : t.get_all_checks_ref()) {
            if (!handle_even_sides(t, c_main)) found_empty = true;
        }
        if (found_empty) continue;
        // Now that is done, we want to add blue checks. This involves finding
        // red-green cycles (blue strings). A blue check will be assigned to
        // the shortest blue strings.
        std::deque<cycle_t> cycles;
        for (const sptr<check_t>& c : t.get_all_checks_ref()) {
            auto arr = get_cycles_on_check(t,c);
            vtils::push_back_range(cycles, arr.begin(), arr.end());
        }
        std::sort(cycles.begin(), cycles.end(),
                [] (auto a, auto b) { return a.size() < b.size(); });
        // Now make the blue checks.
        int k = -1;
        int cycles_applied = 0;
        while (cycles.size()) {
            cycle_t cyc = cycles.front();
            cycles.pop_front();
            // Make sure the cycle is free.
            bool is_free = true;
            std::vector<int> sides;
            for (int s = 0; s < cyc.size(); s++) {
                sptr<check_t> c_curr = cyc.at(s),
                              c_prev = (s==0) ? cyc.back() : cyc.at(s-1),
                              c_next = (s==cyc.size()-1) ? cyc.front() : cyc.at(s+1);
                // Find the appropriate side for c_curr
                int _s = -1;
                for (int i = 1; i < 8; i += 2) {
                    sptr<check_t> x = c_curr->get(i-1),
                                  y = c_curr->get(i+1);
                    if (x == c_prev && y == c_next) {
                        _s = i;
                        break;
                    }
                }
                if (_s < 0 || c_curr->get(_s) != nullptr) {
                    is_free = false;
                    break;
                }
                sides.push_back(_s);
            }
            if (!is_free) {
                continue;
            }

            // Make a blue check.
            sptr<check_t> c_blue = t.add_check_at(2, cyc.size(), k, k);
            for (int s = 0; s < cyc.size(); s++) {
                sptr<check_t> c_curr = cyc.at(s);
                int _s = sides.at(s);
                link(c_blue, s+cyc.front()->color, c_curr, _s);
            }
            cycles_applied++;
            k--;
        }
        if (is_sample_good(t)) {
            samples.push_back(t);
        }
    }
    return samples;
}

bool
MonteCarloManager::is_sample_good(Tiling& t) {
    // Check that every RG side has at least one blue neighbor.
    for (int c = 0; c < 2; c++) {
        for (sptr<check_t> c : t.get_checks_ref(c)) {
            for (int s = 0; s < 8; s += 2) {
                sptr<check_t> left = c->get(s-1),
                              right = c->get(s+1);
                if (left == nullptr && right == nullptr) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool
MonteCarloManager::handle_even_sides(Tiling& t, sptr<check_t>& c_main) {
    for (int side = 0; side < 8; side += 2) {
        if (c_main->get(side) != nullptr) continue;
        // Find candidate with free side.
//      int _side = (side + 2*(rng()%3)+2) % 8;
        int _side = (side+4) % 8;
        auto arr = get_candidates_where(t, c_main, _side, 1-c_main->color);
        if (arr.empty()) return false;  // Bad sample.
        sptr<check_t> c = arr[rng() % arr.size()];
        link(c_main, side, c, _side);
    }
    return true;
}

std::array<cycle_t, 4>
MonteCarloManager::get_cycles_on_check(Tiling& t, sptr<check_t> c_main) {
    std::array<cycle_t, 4> cycles;

    int k = 0;
    for (int s = 0; s < 8; s += 4) {
        sptr<check_t> c_root = c_main->get(s),
                      c_t1 = c_main->get(s-2),
                      c_t2 = c_main->get(s+2);
        std::deque<sptr<check_t>> fifo{ c_root };
        std::set<sptr<check_t>> visited{ c_main, c_root };
        std::map<sptr<check_t>, std::vector<sptr<check_t>>> prev_map;
        while (fifo.size()) {
            sptr<check_t> c = fifo.front();
            fifo.pop_front();

            std::vector<sptr<check_t>> p(prev_map[c]);
            p.push_back(c);
            for (int i = 0; i < 8; i += 2) {
                sptr<check_t> _c = c->get(i);
                if (_c == nullptr) {
                    std::cout << print_check(c) << "." << i << " is null " << std::endl;
                    exit(1);
                }
                if (!visited.count(_c)) {
                    prev_map[_c] = p;
                    fifo.push_back(_c);
                    visited.insert(_c);
                }
            }
        }
        // c_t1 and c_t2 should have paths in prev.
        for (sptr<check_t> c : {c_t1, c_t2}) {
            cycle_t& cyc = cycles[k++];
            cyc.push_back(c_main);
            vtils::push_back_range(cyc, prev_map.at(c));
            cyc.push_back(c);
        }
    }
    return cycles;
}

std::vector<sptr<check_t>>
MonteCarloManager::get_candidates_where(
        Tiling& t, sptr<check_t> root, int side, int color) 
{
    const int LR2 = config.link_radius*config.link_radius;

    const bool root_is_row_edge = (root->i == 0) || (root->i == config.r-1);
    const bool root_is_col_edge = (root->j == 0) || (root->j == config.c-1);

    std::vector<sptr<check_t>> arr( t.get_checks_ref(color) );
    for (auto it = arr.begin(); it != arr.end(); ) {
        const bool is_row_edge = ((*it)->i == 0) || ((*it)->i == config.r-1);
        const bool is_col_edge = ((*it)->j == 0) || ((*it)->j == config.c-1);
        // Check distance from root.
        int di = root->i - (*it)->i,
            dj = root->j - (*it)->j;
        if (root_is_row_edge && is_row_edge) di = 0;
        if (root_is_col_edge && is_col_edge) dj = 0;
        const int r2 = di*di + dj*dj;
        if (r2 > LR2 
//              || (r2 == 1 && (!root_is_row_edge || !root_is_col_edge))
                || (*it)->get(side) != nullptr
                || *it == root
                || (*it)->get(side-2) == root
                || (*it)->get(side+2) == root) 
        {
            it = arr.erase(it);
        } else {
            it++;
        }
    }
    return arr;
}

}   // cgen