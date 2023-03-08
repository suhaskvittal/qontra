/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_h
#define FLEECE_h

#include <stim.h>

#include <PerfectMatching.h>

#include "decoder.h"
#include "defs.h"
#include "fleece/lattice_graph.h"
#include "graph/dijkstra.h"

#include <fstream>
#include <iostream>
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

//#define FLEECE_DEBUG
//#define FLEECE_DEBUG2

namespace qrc {

class Fleece {
public:
    Fleece(const stim::CircuitGenParameters&,
            std::mt19937_64& rng,
            char reset_basis='Z',
            char output_basis='Z');
    ~Fleece();

    stim::simd_bit_table create_syndromes(uint64_t shots, uint disable_leakage_at_round);

    bool fake_run;
private:
    void apply_reset(uint32_t qubit, bool add_error=true);
    void apply_round_start_error(uint32_t qubit);
    void apply_H(uint32_t qubit);
    void apply_CX(uint32_t, uint32_t);
    void apply_measure(uint32_t qubit, bool add_error=true);
    void apply_SWAP(uint32_t, uint32_t, bool add_error=true);

    const stim::CircuitGenParameters circuit_params;
    const char reset_basis;
    const char output_basis;

    stim::Circuit base_circuit;
    stim::FrameSimulator * sim;

#ifdef FLEECE_DEBUG
    std::string op_buffer;
#endif

    fleece::LatticeGraph lattice_graph;
    std::vector<fleece::LatticeGraph::Vertex*> data_qubits;
    std::vector<fleece::LatticeGraph::Vertex*> parity_qubits;

    std::map<uint32_t, fp_t> rounddp_table;
    std::map<uint32_t, fp_t> clifforddp1_table;
    std::map<std::vector<uint32_t>, fp_t> clifforddp2_table;
    std::map<uint32_t, fp_t> premeasflip_table;
    std::map<uint32_t, fp_t> postresflip_table;

    std::map<uint32_t, fp_t> roundleak_table;
    std::map<std::vector<uint32_t>, fp_t> cliffordleak_table;
    std::map<uint32_t, fp_t> postresleak_table;

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
