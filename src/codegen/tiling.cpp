/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

#include "codegen/tiling.h"

#include <iostream>

namespace cgen {

check_t::check_t(int color, int s, int i, int j)
    :color(color),
    sides(s, nullptr),
    i(i),
    j(j)
{}

void
link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s) {
    if (c1->get(c1s) != nullptr || c2->get(c2s) != nullptr) {
        std::cerr << "Writing to already allocated side: "
            << print_check(c1) << "." << c1s << " = " << print_check(c1->get(c1s))
            << ", " << print_check(c2) << "." << c2s << " = "
            << print_check(c2->get(c2s)) << ".\n";
        exit(1);
    }
    c1->get(c1s) = c2;
    c2->get(c2s) = c1;
    c1->get_side_of(c2, c2s) = c1->san(c1s);
    c2->get_side_of(c1, c1s) = c2->san(c2s);
}

Tiling::Tiling(int r, int c)
    :r(r),
    c(c),
    in_plane( r, std::vector<sptr<check_t>>(c, nullptr) ),
    checks_by_color(),
    all()
{
    for (int j = 0; j < c; j++) {
        for (int i = 0; i < r; i++) {
            sptr<check_t> c1 = add_check_at((i+j)%2, 8, i, j);
            if (i > 0) {
                // Create link to previous check.
                sptr<check_t> c2 = at(i-1, j);
                link(c1, 0, c2, 4);
            }
        }
        if (j == 0) continue;
        for (int i : {0,r-1}) {
            // Create link to check in previous column across row i.
            sptr<check_t> c1 = at(i,j),
                          c2 = at(i,j-1);
            link(c1, 6, c2, 2);
        }
    }
}

sptr<check_t>
Tiling::add_check_at(int color, int sides, int i, int j) {
    sptr<check_t> c = std::make_shared<check_t>(color, sides, i, j);
    if (i >= 0 || j >= 0) {
        in_plane[i][j] = c;
    }
    checks_by_color[color].push_back(c);
    all.push_back(c);
    return c;
}

}   // cgen
