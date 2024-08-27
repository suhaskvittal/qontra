/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#ifndef CGEN_SC_MC_h
#define CGEN_SC_MC_h

#include "codegen/surface_codes/utils.h"
#include "codegen/surface_codes/coord.h"

#include <vtils/timer.h>

#include <array>
#include <map>
#include <random>
#include <tuple>
#include <vector>
#include <utility>

namespace cgen {

template <int R, int S, size_t W, size_t H>
struct StarClass {
    static constexpr size_t TOTAL_COUNT = W*H;

    StarClass(Star<R,S>* r, coord_t rloc)
        :repr(r),
        repr_loc(rloc),
        coord_map(),
        elements()
    {
        // Setup coord_map.
        elements[0] = r;
        coord_map[r] = PZERO;
        for (uint8_t i = 0; i < static_cast<uint8_t>(W); i++) {
            for (uint8_t j = 0; j < static_cast<uint8_t>(H); j++) {
                if (i == 0 && j == 0) continue;
                Star<R,S>* check = new Star<R,S>;
                coord_map[check] = make_coord(i,j);
                elements[i*H+j] = check;
            }
        }
    }

    void connect_with_delta(StarClass<R,S,W,H>*, const coord_t& delta, uint8_t qi, uint8_t qj);

    Star<R,S>* repr;
    coord_t repr_loc;
    std::map<Star<R,S>*, coord_t> coord_map;
    std::array<Star<R,S>*, TOTAL_COUNT> elements;
};

template <int R, int S, 
            size_t GRID_WIDTH,
            size_t GRID_HEIGHT,
            size_t CELL_WIDTH=2, 
            size_t CELL_HEIGHT=2>
class MonteCarloManager {
public:
    static constexpr size_t CBW = GRID_WIDTH/CELL_WIDTH;
    static constexpr size_t CBH = GRID_HEIGHT/CELL_HEIGHT;

    MonteCarloManager()
        :fpdst(0.0, 1.0),
        rng(0)
    {}

    typedef StarClass<R,S,CBW,CBH>  S_eqc_t;

    struct sample_t {
        ~sample_t() {
            for (S_eqc_t* x : star_classes) delete x;
            for (Star<R,S>* x : star_checks) delete x;
        }

        coord_t get_loc(Star<R,S>* x) {
            S_eqc_t* c_p = s_eqc_map.at(x);
            return c_add(c_p->repr_loc, c_p->coord_map.at(x));
        }

        std::map< Star<R,S>*, S_eqc_t* > s_eqc_map;
        // Holding structures:
        std::vector< S_eqc_t* > star_classes;
        std::vector< Star<R,S>* > star_checks;
    };

    sample_t make_sample(void);

    struct {
        fp_t star_init_prob = 0.7;

        size_t max_link_len = 12;
    } config;

    struct {
        uint64_t total_time_in_init = 0.0;
        uint64_t total_time_overall = 0.0;
    } stats;
private:
    sample_t init_sample(void);

    std::vector< std::tuple<Star<R,S>*, S_eqc_t*, coord_t> >  // candidate, class, delta
        get_candidates_for_star_check(Star<R,S>* rt, const coord_t& rt_loc, const sample_t&);

    std::uniform_real_distribution<> fpdst;
    std::mt19937_64 rng;

    vtils::Timer global_timer;
    vtils::Timer local_timer;
};

}   // cgen

#include "mc.inl"

#endif  // CGEN_SC_MC_h
