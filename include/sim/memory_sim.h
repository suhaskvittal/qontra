/*
 *  author: Suhas Vittal
 *  date:   10 December 2023
 * */

#ifndef SIM_SURFACE_CODE_MEMORY_h
#define SIM_SURFACE_CODE_MEMORY_h

#include "experiments.h"
#include "instruction.h"
#include "graph/lattice_graph.h"
#include "sim/frame_sim.h"
#include "tables.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <stim.h>

#include <math.h>
#include <mpi.h>
#include <stdio.h>

namespace qontra {
namespace sim {

graph::LatticeGraph surface_code_lattice_graph(uint distance);

}   // sim

class MemorySimulator {
public:
    MemorySimulator(graph::LatticeGraph&);
    ~MemorySimulator();

    void    reset();

    void    run(uint64_t shots);

    struct {
        uint distance;
        uint rounds;
        ErrorTable  errors;
        TimeTable   timing;

        std::string stim_output_file;
        std::string syndrome_output_folder;
        std::string data_output_file;

        bool is_memory_x=false;
    } config;
private:
    typedef fp_t time_t;

    void    run_batch(uint64_t shots_in_batch);

    time_t  do_gate(std::string op, std::vector<uint> operands);
    time_t  do_measurement(std::vector<uint> operands);
    void    create_event_or_obs(std::vector<uint> operands, bool create_event=true);
    void    inject_timing_error(std::vector<uint> qubits);
    void    inject_idling_error_positive(std::vector<uint> on_qubits);
    void    inject_idling_error_negative(std::vector<uint> not_on_qubits);
    
    graph::LatticeGraph lattice_graph;
    uint64_t n_qubits;
    uint64_t n_detection_events;
    uint64_t n_observables;

    std::vector<uint>   all_qubits;
    std::vector<uint>   data_qubits;
    std::vector<uint>   parity_qubits;
    std::vector<uint>   xp_qubits;
    std::vector<uint>   zp_qubits;

    std::vector<uint>   syndrome_parity_qubits;

    std::vector<std::vector<uint>>  obs_list;

    FrameSimulator* sim;
    stim::simd_bit_table syndromes;
    // Tracking data structures.
    time_t  elapsed_time;
    std::map<uint64_t, time_t> shot_time_delta_map;

    uint64_t meas_ctr;
    uint64_t event_ctr;
    uint64_t obs_ctr;
    std::map<uint, uint64_t> meas_ctr_map;
    std::map<uint, uint64_t> event_ctr_map;

    stim::Circuit   sample_circuit;
    bool            is_recording_stim_instructions;
};

}   // qontra

#endif  // SIM_SURFACE_CODE_MEMORY_h
