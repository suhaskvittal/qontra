/*
 *  author: Suhas Vittal
 *  date:   20 June 2023
 * */

#ifndef CONTROL_SIM_h
#define CONTROL_SIM_h

#include "decoder/decoder.h"
#include "defs.h"
#include "experiments.h"
#include "instruction.h"
#include "sim/clifford_sim.h"
#include "tables.h"

#include <functional>
#include <map>
#include <vector>
#include <random>

#include <math.h>

namespace qontra {

// The ControlSimulator simulates high-level aspects
// of the control architecture, mainly:
//  (1) Instruction scheduling of QUANTUM instructions.
//  (2) Decoding and tracking of Pauli frames

class ControlSimulator {
public:
    ControlSimulator(uint n_qubits,
                    decoder::Decoder* dec,
                    const schedule_t& program)
        // Stats
        :latency(experiments::G_SHOTS_PER_BATCH, 64),
        sim_time(0),
        // Simulation state tracking
        decoder_busy(experiments::G_SHOTS_PER_BATCH, 64),
        trial_done(experiments::G_SHOTS_PER_BATCH),
        // Simulation structures
        decoder(dec),
        csim(n_qubits, experiments::G_SHOTS_PER_BATCH),
        pauli_frames(experiments::G_SHOTS_PER_BATCH, 128),
        event_history(4096, experiments::G_SHOTS_PER_BATCH),
        obs_buffer(128, experiments::G_SHOTS_PER_BATCH),
        obs_buffer_max_written(0),
        // Microarchitecture
        program(program),
        pc(experiments::G_SHOTS_PER_BATCH, 64),
        // IF io
        if_stall(experiments::G_SHOTS_PER_BATCH),
        if_pc(experiments::G_SHOTS_PER_BATCH, 64),
        if_id_valid(experiments::G_SHOTS_PER_BATCH),
        // ID io
        id_stall(experiments::G_SHOTS_PER_BATCH),
        id_pc(experiments::G_SHOTS_PER_BATCH, 64),
        id_qex_valid(experiments::G_SHOTS_PER_BATCH),
        // QEX io
        qex_stall(experiments::G_SHOTS_PER_BATCH),
        qex_rt_valid(experiments::G_SHOTS_PER_BATCH),
        qex_pc(experiments::G_SHOTS_PER_BATCH, 64),
        qex_qubit_busy(experiments::G_SHOTS_PER_BATCH, 64*n_qubits),
        // RT io
        rt_stall(experiments::G_SHOTS_PER_BATCH),
        // Other
        n_qubits(n_qubits),
        rng(0)
    {}

    void run(uint64_t shots);
    void clear();

    struct params_t {
        // Simulation parameters
        uint64_t    kill_batch_after_time_elapsed = 1000L * 1'000'000'000L; 
                                            // In nanoseconds, halts
                                            // the simulator after this
                                            // much walltime has elapsed.
        uint        verbose = 0;
        // Control system configuration
        fp_t        clock_frequency = 250e6;
        bool        decoder_is_ideal = true;    // Decoder only takes 1ns to run.
        // Configuration of quantum computer
        TimeTable   timing;
        ErrorTable  errors;

        uint64_t    apply_periodic_errors_at_t = 1000;
                                            // Apply state errors on all qubits
                                            // as frequently as shown.
        bool        simulate_periodic_as_dpo_and_dph = true;
                                            // If false, we just do depolarizing
                                            // errors instead.
    };
    params_t params;
    // Stats: add any new stats if necessary.
    stim::simd_bit_table                latency;
                                        // Latency is a 64-bit integer: 
                                        // units are ns.
    std::map<vlw_t, uint64_t>           prob_histograms;
                                        // Stores the number of shots that 
                                        // give certain output measurements.
    uint64_t    n_trials_killed;
    uint64_t    latency_max;
    uint64_t    latency_mean;
    uint64_t    latency_std;

    uint64_t    sim_time;
private:
    void        IF(void);   // Instruction fetch
    void        ID(void);   // Instruction decode
    void        QEX(void);  // Quantum instruction execution
    void        RT(void);   // Retire
    
    typedef std::function<void(std::vector<uint>, std::vector<fp_t>)>   ef_t;

    void        apply_gate_error(Instruction&);
    void        apply_periodic_error(void);

    uint64_t    shots_in_curr_batch;
    Timer       timer;

    // Simulation state tracking
    stim::simd_bit_table    decoder_busy;
    stim::simd_bits         trial_done; // Active-low: if 1, then not finished.
    // Simulation structures:
    decoder::Decoder*       decoder; // Should return XZ frame changes
    CliffordSimulator       csim;
    stim::simd_bit_table    pauli_frames;
    stim::simd_bit_table    event_history;
    stim::simd_bit_table    obs_buffer;
    uint64_t                obs_buffer_max_written;
    // Microarchitecture
    schedule_t              program;    // Program should remain constant
    stim::simd_bit_table    pc;

    // IF io
    stim::simd_bits         if_stall;
    stim::simd_bits         if_id_valid;
    stim::simd_bit_table    if_pc;
    // ID io
    stim::simd_bits         id_stall;
    stim::simd_bits         id_qex_valid;
    stim::simd_bit_table    id_pc;
    // QEX io
    stim::simd_bits         qex_stall;
    stim::simd_bits         qex_rt_valid;
    stim::simd_bit_table    qex_pc;
    stim::simd_bit_table    qex_qubit_busy;
    // RT io
    stim::simd_bits         rt_stall;

    const uint      n_qubits;
    std::mt19937_64 rng;
};

}   // qontra

#endif  // CONTROL_SIM_h
