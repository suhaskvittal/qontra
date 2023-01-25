/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_h
#define FLEECE_h

#include <stim.h>

#include <PerfectMatching.h>

#include "defs.h"
#include "fleece/lattice_graph.h"
#include "graph/dijkstra.h"

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
            char output_basis='Z');

    typedef std::pair<stim::simd_bit_table, stim::simd_bit_table> SyndromeOutput;

    SyndromeOutput generate_syndromes(
            uint64_t shots, uint last_leakage_round, uint64_t seed=0);
    std::map<uint32_t, uint8_t> get_data_qubit_states(uint64_t shot_number);

    bool fake_run;
private:
    void correct_leak(uint64_t shots);
    void correct_data_qubit(uint64_t shot_number, 
            uint detector, fleece::LatticeGraph::Vertex*);
    void correct_parity_qubit(uint64_t shot_number,
            uint detector, fleece::LatticeGraph::Vertex*);

    uint get_level(uint);
    uint jump_to_level(uint base, uint level);
    uint base_detector(uint);
    uint prev_detector(uint);
    uint next_detector(uint);
    uint32_t get_measurement_time(uint);

    const stim::Circuit circuit;
    const uint detectors_per_round;
    const bool double_stabilizer;
    const char reset_basis;
    const char output_basis;

    stim::FrameSimulator sim;
    stim::simd_bit_table meas_results;
    stim::simd_bit_table leak_results;

    fleece::LatticeGraph lattice_graph;
    graph::PathTable<fleece::LatticeGraph::Vertex> path_table;

    uint curr_min_detector;
    uint curr_max_detector;

    std::set<uint32_t> await_detector_set;

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
