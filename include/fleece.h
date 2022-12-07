/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_h
#define FLEECE_h

#include <stim.h>

#include "defs.h"
#include "fleece/lattice_graph.h"

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <utility>

/*
 *  FLEECE = Fast LEakagE CorrEction
 * */

namespace qrc {

const uint64_t MAX_SHOTS = 100000;

class Fleece {
public:
    Fleece(const stim::Circuit&,
            std::mt19937_64& rng,
            uint detectors_per_round,
            bool double_stabilizer=true,
            char reset_basis='Z',
            char output_basis='Z',
            uint32_t tower_cutoff=2);

    typedef std::pair<stim::simd_bit_table, stim::simd_bit_table> SyndromeOutput;

    SyndromeOutput generate_syndromes(uint64_t shots);
private:
    void tower_correct(uint64_t shots);

    uint base_detector(uint);
    uint next_detector(uint);

    const stim::Circuit circuit;
    const uint detectors_per_round;
    const bool double_stabilizer;
    const char reset_basis;
    const char output_basis;
    const uint32_t tower_cutoff;

    stim::FrameSimulator sim;
    stim::simd_bit_table meas_results;
    stim::simd_bit_table leak_results;

    fleece::LatticeGraph lattice_graph;

    uint curr_min_detector;
    uint curr_max_detector;

    std::mt19937_64 rng;
};

} // qrc

#endif
