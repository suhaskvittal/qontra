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
#include "sim/components/syndrome_predictor.h"
#include "sim/frame_sim.h"
#include "tables.h"

#include <algorithm>
#include <filesystem>
#include <functional>
#include <map>
#include <vector>
#include <random>

#include <math.h>
#include <stdio.h>

namespace qontra {

// The ControlSimulator simulates high-level aspects
// of the control architecture, mainly:
//  (1) Instruction scheduling of QUANTUM instructions.
//  (2) Decoding and tracking of Pauli frames

namespace ctrlsim {

extern uint64_t     EVENT_HISTORY_SIZE; // Default is 1024*8 (1KB)
extern uint64_t     OBS_BUFFER_SIZE;    // Default is 128.
extern uint64_t     TRACES_PER_SHOT;    // Default is 1.

}

class ControlSimulator {
public:
    ControlSimulator(uint n_qubits, const schedule_t&, StateSimulator*);

    void run(uint64_t shots);
    void clear();

    void            build_error_model();
    stim::Circuit   get_error_model(void) { return canonical_circuit; }

    void    load_simulator(StateSimulator* ss) { qsim = ss; }
    void    load_decoder(decoder::Decoder* dec) { decoder = dec; }
    void    load_predictor(sim::SyndromePredictor* pr) 
                { syndrome_predictor = pr; }

    struct {
        // Simulation parameters
        uint64_t    kill_batch_after_time_elapsed = 1000; 
                                            // In seconds, halts
                                            // the simulator after this
                                            // much walltime has elapsed.
        uint        verbose = 0;
        // Fast-forwarding parameters
        //      Fast-forwarding the simulator makes all instructions take
        //      1 cycle and have no execution latency (i.e. a measurement
        //      now takes 0ns). This is best used for quick simulations.
        //
        //      Fast-forwarding can be integrated into the simulation via
        //      virtual instructions ffstart and ffend (see instruction.h).
        bool        enable_fast_forwarding = true;
        uint64_t    ff_apply_periodic_errors_at_t = 250;
        uint64_t    ff_periodic_error_assume_time_elapsed = 1000;
        // Control system configuration
        fp_t        clock_frequency = 250e6;
        bool        decoder_is_ideal = true;    // Decoder only takes 1ns to run.

        // Configuration of quantum computer
        TimeTable   timing;
        ErrorTable  errors;
        uint64_t    apply_periodic_errors_at_t = 100;
                                            // Apply state errors on all qubits
                                            // as frequently as shown.
        bool        simulate_periodic_as_dpo_and_dph = false;
                                            // If false, we just do depolarizing
                                            // errors instead.
        // Saving syndrome data
        bool                save_syndromes_to_file = false;
        std::string         syndrome_output_folder; 
                                        // Syndromes are saved in Stim's
                                        // .dets file format. One file is
                                        // created for each batch.
        std::vector<uint>   save_observables{0};
                                        // A list of observable indices to
                                        // save.
    } params;
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
    typedef std::array<Instruction, 4096>   imem;

    void        IF(void);   // Instruction fetch
    void        ID(void);   // Instruction decode
    void        QEX(void);  // Quantum instruction execution
    void        RT(void);   // Retire
    
    typedef std::function<void(std::vector<uint>, std::vector<fp_t>)>   ef_t;

    void        apply_gate_error(Instruction&);
    void        apply_periodic_error(fp_t t);

    uint64_t    shots_in_curr_batch;
    Timer       timer;

    // Simulation state tracking
    stim::simd_bit_table    decoder_busy;
    stim::simd_bits         trial_done; // Active-low: if 1, then not finished.
    // Simulation structures:
    decoder::Decoder*       decoder;
    StateSimulator*         qsim;
    stim::simd_bit_table    pauli_frames;
    stim::simd_bit_table    event_history;
    stim::simd_bit_table    obs_buffer;
    uint64_t                obs_buffer_max_written;

    stim::simd_bit_table    event_trace;    // Stores events to be written
                                            // if specified in the params.
    stim::simd_bit_table    event_trace_index;
    uint64_t                event_trace_max_written;
    // Microarchitecture
    imem                    program;    // Program should remain constant
    stim::simd_bit_table    pc;

    //
    //  Data structures for Selective Syndrome Extraction
    //
    stim::simd_bit_table    sig_m_spec;     // Active if the measurement was
                                            // speculated.
    stim::simd_bit_table    val_m_spec;     // Value of speculated measurement.

    sim::SyndromePredictor* syndrome_predictor;

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
    stim::simd_bits         qex_br_taken;
    // RT io
    stim::simd_bits         rt_stall;

    const uint      n_qubits;
    std::mt19937_64 rng;

    // Canonical circuit structures.
    bool                    is_building_canonical_circuit;
    stim::Circuit           canonical_circuit;
    std::map<uint, uint>    measurement_order_to_loc;
    uint                    measurement_count;
    
    // Fast forwarding structures.
    bool    is_fast_forwarding;
    bool    apply_pending_errors;
};

}   // qontra

#endif  // CONTROL_SIM_h
