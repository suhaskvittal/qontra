/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

#ifndef CODEGEN_TILING_h
#define CODEGEN_TILING_h

#include <qontra/defs.h>

#include <vtils/mat2.h>
#include <vtils/two_level_map.h>

#include <array>
#include <string>
#include <vector>

// CONSISTENCY RULES:
//  Red checks --> even sides = green, odd sides = blue
//  Green checks --> even sides = red, odd sides = green
//  Blue checks --> even sides = red, odd sides = green

namespace cgen {

struct check_t {
    check_t(int color, int sides, int i, int j);

    inline sptr<check_t>& get(int s) {
        return sides[san(s)];
    }

    inline int& get_side_of(sptr<check_t> c, int cs) {
        return side_map[c][c->san(cs)];
    }

    inline const int& get_const_side_of(sptr<check_t> c, int cs) const {
        return side_map.at(c).at(c->san(cs));
    }

    inline int san(int s) const {
        if (s < 0) return size() + s;
        else       return s % size();
    }

    inline int size() const { return sides.size(); }

    int color;
    int i;
    int j;
    std::vector<sptr<check_t>> sides;
    vtils::TwoLevelMap<sptr<check_t>, int, int> side_map;
};

inline std::string print_check(sptr<check_t> c) {
    if (c == nullptr) return "NIL";
    else {
        std::string color;
        if (c->color == 0) color = "r";
        else if (c->color == 1) color = "g";
        else color = "b";
        return color + "(" + std::to_string(c->i) + "," + std::to_string(c->j) + ")";
    }
}

void link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s);

class Tiling {
public:
    Tiling(int r, int c);
    Tiling(const Tiling&) =default;

    sptr<check_t>& operator()(int i, int j);
    sptr<check_t> at(int i, int j) const;

    sptr<check_t> add_check_at(int color, int sides, int i=-1, int j=-1);

    std::vector<sptr<check_t>> get_checks(int color) const;
    std::vector<sptr<check_t>> get_all_checks(void) const;

    const std::vector<sptr<check_t>>& get_checks_ref(int color) const;
    const std::vector<sptr<check_t>>& get_all_checks_ref(void) const;

    const int r;
    const int c;
private:
    std::vector<std::vector<sptr<check_t>>> in_plane;
    std::array<std::vector<sptr<check_t>>, 3> checks_by_color;
    std::vector<sptr<check_t>> all;
};

}   // cgen

#include "inl/tiling.inl"

#endif  // CODEGEN_TILING_h
