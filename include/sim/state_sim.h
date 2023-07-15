/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 *
 *  A virtual class for a State Simulator.
 * */

#ifndef STATE_SIM_h
#define STATE_SIM_h

#include "defs.h"

#include <random>
#include <vector>

#include <stim.h>

namespace qontra {

namespace statesim {

extern uint64_t G_RECORD_SPACE_SIZE;

}   //  statesim

// This class is generic, and intended to be extensible.
//
// Right now, QontraSim implements the following simulators:
//  (1) FrameSimulator, essentially a clone of Stim's FrameSimulator
//      but more adopted for our purposes. Use this for memory experiments
//      or other experiments where we just care about Pauli frame
//      changes.
//  (2) CliffordSimulator, an optimized simulator. Use this when simulating
//      proper quantum programs.

class StateSimulator {
public:
    StateSimulator(uint n, uint64_t max_shots)
        :n_qubits(n),
        max_shots(max_shots),
        record_table(statesim::G_RECORD_SPACE_SIZE, max_shots),
        lock_table(n, max_shots),
        record_table_cpy(statesim::G_RECORD_SPACE_SIZE, max_shots),
        lock_table_cpy(n, max_shots),
        rng(0)
    {
        reset_sim();
    }

    void    set_seed(uint64_t x) { rng.seed(x); }

    virtual void    reset_sim(void) {
        record_table.clear();
        lock_table.clear();
    }

    virtual void    H(std::vector<uint>) =0;
    virtual void    X(std::vector<uint>) =0;
    virtual void    Z(std::vector<uint>) =0;
    virtual void    S(std::vector<uint>) =0;
    virtual void    CX(std::vector<uint>) =0;
    virtual void    M(std::vector<uint>, int record=-1) =0;
    virtual void    R(std::vector<uint>) =0;

    virtual void    eDPO(std::vector<uint>, std::vector<fp_t>) =0;
    virtual void    eDPH(std::vector<uint>, std::vector<fp_t>) =0;

    virtual void    eDP1(std::vector<uint>, std::vector<fp_t>) =0;
    virtual void    eDP2(std::vector<uint>, std::vector<fp_t>) =0;
    virtual void    eX(std::vector<uint>, std::vector<fp_t>) =0;
    virtual void    eLI(std::vector<uint>, std::vector<fp_t>) =0;
    virtual void    eLT(std::vector<uint>, std::vector<fp_t>) =0;

    void    shift_record_by(uint64_t);

    virtual void    snapshot(void);
                            // Saves the current state of the simulator.
    virtual void    rollback_where(stim::simd_bits_range_ref);
                            // Rolls back the state to the snapshot

    stim::simd_bit_table    record_table;
    stim::simd_bit_table    lock_table;
    uint64_t                shots;
protected:
    stim::simd_bit_table    record_table_cpy;
    stim::simd_bit_table    lock_table_cpy;

    std::mt19937_64 rng;

    const uint      n_qubits;
    const uint64_t  max_shots;
};

void    copy_where(stim::simd_bits_range_ref from,
                    stim::simd_bits_range_ref to,
                    stim::simd_bits_range_ref pred);

}   // qontra

#endif  // STATE_SIM_h
