/*
 *  author: Suhas Vittal
 *  date:   8 August 2024
 * */

#ifndef CODEGEN_MC_h
#define CODEGEN_MC_h

#include "codegen/tiling.h"

#include <map>
#include <random>

namespace cgen {

#define GET_RAND_IT(arr, rng) arr.begin() + rng()

typedef std::vector<sptr<check_t>> cycle_t;

class MonteCarloManager {
public:
    MonteCarloManager(int seed);

    std::vector<Tiling> run(uint64_t samples);

    struct {
        int r = 4;
        int c = 4;
    } config;
private:
    bool is_sample_good(Tiling&);

    bool handle_even_sides(Tiling&, sptr<check_t>&);
    std::array<cycle_t, 4> get_cycles_on_check(Tiling&, sptr<check_t>);

    // Gets all candidates
    std::vector<sptr<check_t>>
        get_candidates_where(Tiling&, sptr<check_t> root, int side, int color);
    
    std::mt19937_64 rng;
};

}   // cgen

#include "inl/monte_carlo.inl"

#endif  // CODEGEN_MC_h
