/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_h
#define FLEECE_h

#include <stim.h>

#include "defs.h"

#include <map>
#include <random>
#include <string>
#include <vector>


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
            char output_basis='Z');

    stim::simd_bit_table generate_syndromes(uint64_t shots);
private:
    void tower_correct(uint64_t shots);

    uint base_detector(uint);

    const stim::Circuit circuit;
    const uint detectors_per_round;
    const bool double_stabilizer;
    const char reset_basis;
    const char output_basis;

    stim::FrameSimulator sim;
    stim::simd_bit_table current_results;

    DecodingGraph decoding_graph;
    PathTable path_table;

    std::map<uint, std::vector<uint>> stabilizer_to_data;

    uint curr_min_detector;
};

} // qrc

#endif
