/* author: Suhas Vittal
 *  date:   7 August 2024
 * */

#ifndef CODEGEN_TILING_h
#define CODEGEN_TILING_h

#include <qontra/defs.h>

#include <vtils/mat2.h>
#include <vtils/two_level_map.h>

#include <array>
#include <vector>
#include <iostream>

// CONSISTENCY RULES:
//  Red checks --> even sides = green, odd sides = blue
//  Green checks --> even sides = red, odd sides = green
//  Blue checks --> even sides = red, odd sides = green

namespace cgen {

struct check_t {
    check_t(int color, int sides, int i, int j);

    inline sptr<check_t>& get(int s) {
        if (s < 0) return sides[sides.size()-s];
        else       return sides[s % sides.size()];
    }

    inline int size() const { return sides.size(); }

    int color;
    int i;
    int j;
    std::vector<sptr<check_t>> sides;
    vtils::TwoLevelMap<sptr<check_t>, int, int> side_map;
};

void link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s);

class Tiling {
public:
    Tiling(int r, int c, int blue_sides);
    Tiling(const Tiling&) =default;

    std::vector<Tiling> generate_sample_tilings(uint64_t samples, int seed) const;

    vtils::Mat2 to_matrix(void) const;

    sptr<check_t>& operator()(int i, int j);
    sptr<check_t> at(int i, int j) const;

    sptr<check_t> add_check_at(int color, int sides, int i=-1, int j=-1);
    std::vector<sptr<check_t>> get_all_checks(void) const;

    const int r;
    const int c;
    const int blue_sides;
private:
    std::vector<std::vector<sptr<check_t>>> in_plane;
    std::vector<sptr<check_t>> out_of_plane;

    std::vector<sptr<check_t>> all;
};

}   // cgen

#include "inl/tiling.inl"

#endif  // CODEGEN_TILING_h
