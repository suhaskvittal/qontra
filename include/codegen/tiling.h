/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

#ifndef CODEGEN_TILING_h
#define CODEGEN_TILING_h

#include <qontra/defs.h>

#include <array>
#include <vector>

// CONSISTENCY RULES:
//  Red checks --> even sides = green, odd sides = blue
//  Green checks --> even sides = red, odd sides = green
//  Blue checks --> even sides = red, odd sides = green

namespace cgen {

struct check_t {
    check_t(int color);

    int color;
    std::array<sptr<check_t>, 8> sides;
};

void link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s);

class Tiling {
public:
    Tiling(int, int);
    Tiling(const Tiling&) =default;

    std::vector<Tiling> generate_sample_tilings(uint64_t samples, int seed) const;

    sptr<check_t>& operator()(int i, int j);
    sptr<check_t> at(int i, int j) const;

    sptr<check_t> add_check_at(int color, int i=-1, int j=-1);
    std::vector<check_t> get_all_checks(void) const;

    const int r;
    const int c;
private:
    std::vector<std::vector<sptr<check_t>>> in_plane;
    std::vector<sptr<check_t>> out_of_plane;

    std::vector<sptr<check_t>> all;
};

std::vector<Tiling>
generate_sample_tilings_from_base(uint64_t samples, const Tiling&, int seed);

}   // cgen

#endif  // CODEGEN_TILING_h
