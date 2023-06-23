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

namespace qontra {

// The ControlSimulator simulates high-level aspects
// of the control architecture, mainly:
//  (1) Instruction scheduling of QUANTUM instructions.
//  (2) Decoding and tracking of Pauli frames
//

class ControlSimulator {
    ControlSimulator(uint n_qubits,
                    decoder::Decoder* dec,
                    const schedule_t<cp::Instruction>& program)
        // Simulation state tracking
        :decoder_busy(experiments::G_SHOTS_PER_BATCH, 64),
        time(0),
        // Simulation structures
        decoder(dec)
        csim(n_qubits, experiments::G_SHOTS_PER_BATCH),
        event_history(8192, experiments::G_SHOTS_PER_BATCH),
        pauli_frames(128, experiments::G_SHOTS_PER_BATCH),
        // Microarchitecture
        program(program),
        pc(experiments::G_SHOTS_PER_SHOTS, 64),
        // IF io
        if_stall(experiments::G_SHOTS_PER_BATCH),
        if_pc(experiments::G_SHOTS_PER_BATCH, 64),

    {}

    void run(uint64_t shots);
    void clear();

    struct params_t {
        fp_t clock_frequency    = 250e6;
    };
    // Stats: add any new stats if necessary.
    stim::simd_bit_table            latency(experiments::G_SHOTS_PER_BATCH, 64);
                                        // Latency is a 64-bit integer: 
                                        // units are ns.
    std::map<uint64_t, uint64_t>    prob_histograms;
                                        // Stores the number of shots that 
                                        // give certain output measurements.
private:
    void        IF(void);   // Instruction fetch
    void        ID(void);   // Instruction decode
    void        QEX(void);  // Quantum instruction execution
    void        AM(void);   // Adjust measurement

    uint64_t    shots_in_curr_batch;

    // Simulation state tracking
    stim::simd_bit_table    decoder_busy;
    fp_t                    time;
    // Simulation structures:
    Decoder*                decoder; // Should return XZ frame changes
    CliffordSimulator       csim;
    stim::simd_bit_table    event_history;
    stim::simd_bit_table    pauli_frames;   // Two bits per logical qubit.
    // Microarchitecture
    schedule_t<cp::Instruction> program;    // We need not worry about
                                            // branches.
    stim::simd_bit_table        pc;

    // IF io
    stim::simd_bits         if_stall;
    stim::simd_bit_table    if_pc;
    stim::simd_bits         if_id_valid;
    // ID io
    stim::simd_bit_table    id_stall;
    stim::simd_bit_table    id_r1;
    stim::simd_bit_table    id_r2;
    stim::simd_bit_table    id_rd;
    stim::simd_bit_table    id_immd;
};

}   // qontra

#endif  // CONTROL_SIM_h
