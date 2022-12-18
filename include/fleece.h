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
 *
 *  TODO: Implement basis-aware FLEECE. Only supports
 *  Z basis operations at the moment.
 * */

namespace qrc {

class Fleece {
public:
    Fleece(const stim::Circuit&,
            std::mt19937_64& rng,
            uint detectors_per_round,
            bool double_stabilizer=true,
            char reset_basis='Z',
            char output_basis='Z',
            uint32_t tower_cutoff=2,
            uint32_t sliding_window_size=4);

    typedef std::pair<stim::simd_bit_table, stim::simd_bit_table> SyndromeOutput;

    SyndromeOutput generate_syndromes(uint64_t shots, uint64_t seed=0);

    bool toggle_fake_run(void);
private:
    void tower_correct(uint64_t shots);
    void clean_parity_qubit(uint32_t, uint lower_level, uint upper_level, uint64_t shot);

    uint base_detector(uint);
    uint next_detector(uint);

    const stim::Circuit circuit;
    const uint detectors_per_round;
    const bool double_stabilizer;
    const char reset_basis;
    const char output_basis;
    const uint32_t tower_cutoff;
    const uint32_t sliding_window_size;

    stim::FrameSimulator sim;
    stim::simd_bit_table meas_results;
    stim::simd_bit_table leak_results;

    fleece::LatticeGraph lattice_graph;

    uint curr_min_detector;
    uint curr_max_detector;

    bool fake_run;

    std::mt19937_64 rng;
};

stim::simd_bit_table
double_stabilizer_to_single_stabilizer(
        stim::simd_bit_table,
        uint code_dist,
        uint num_detectors,
        uint num_observables,
        uint num_shots,
        bool is_memory_z=true);

} // qrc

#endif
