/*
 *  author: Suhas Vittal
 *  date:   17 August 2023
 * */

#ifndef SIM_MANAGER_h
#define SIM_MANAGER_h

#include "instruction.h"
#include "sim/state_sim.h"
#include "experiments.h"

#include <iostream>

namespace qontra {

typedef std::map<vlw_t, uint64_t>   histogram_t;

inline void
print_histogram(std::ostream& os, uint64_t n, const histogram_t& h, fp_t filter=0.0) {
    // First, get total number of shots.
    uint64_t total_shots = 0;
    for (const auto pair : h) {
        total_shots += pair.second;
    }
    for (const auto pair : h) {
        vlw_t x = pair.first;
        uint64_t c = pair.second;

        fp_t pr = ((fp_t)c) / ((fp_t)total_shots);
        if (pr <= filter)   continue;

        uint64_t k = 0;
        for (uint64_t bits : x) {
            for (uint j = 0; j < 64; j++) {
                if (k == n)  goto print_histogram_done_with_entry;
                std::cout << ((bits >> j) & 1);
                k++;
            }
        }
print_histogram_done_with_entry:
        std::cout << "\t" << pr << "\n";
    }
}

const uint EVENT_HISTORY_SIZE = 4096;
const uint OBS_BUFFER_SIZE = 4096;
const uint PAULI_FRAME_SIZE = 4096;

class SimManager {
public:
    SimManager(const schedule_t& prog)
        :program(prog),
        n_qubits(get_number_of_qubits(prog)),
        // Default initializations:
        decoder(nullptr),
        sim(nullptr),
        timing_table(),
        event_history_table(EVENT_HISTORY_SIZE, experiments::G_SHOTS_PER_BATCH),
        obs_buffer_table(OBS_BUFFER_SIZE, experiments::G_SHOTS_PER_BATCH),
        pauli_frame_table(PAULI_FRAME_SIZE, experiments::G_SHOTS_PER_BATCH),
        meas_record_ctr(0),
        max_obs_written(0)
    {}

    // This is just a wrapper for the corresponding call in
    // experiments.h, which is much more flexible.
    struct result_t {
        histogram_t probability_histogram;
        uint64_t obscnt;
    };

    result_t evaluate_monte_carlo(uint64_t shots);

    struct {
        // Simulation Parameters
        ErrorTable  errors; // Error basic_error_rates for each gate.
        TimeTable   timing; // Timing information

        bool always_inject_timing_errors = false;
        bool ignore_all_errors = false;
    } params;

    Decoder*        decoder;
    StateSimulator* sim;
private:
    void    simulate_batch(uint64_t shots);

    fp_t    get_operation_latency(Instruction);
    void    inject_operation_error(Instruction);
    void    inject_timing_error(fp_t);

    // Structures for simulation:
    std::map<uint64_t, fp_t> timing_table;

    stim::simd_bit_table event_history_table;
    stim::simd_bit_table obs_buffer_table;
    stim::simd_bit_table pauli_frame_table;

    uint64_t meas_record_ctr;
    uint64_t max_obs_written;

    // Program to simulate.
    const schedule_t    program;
    const uint          n_qubits;
};

}   // qontra

#endif  // SIM_MANANGER_h
