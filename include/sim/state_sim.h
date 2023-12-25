/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 *
 *  A virtual class for a State Simulator.
 * */

#ifndef STATE_SIM_h
#define STATE_SIM_h

#include "defs.h"
#include "stimext.h"

#include <random>
#include <tuple>
#include <vector>

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

template <class T, class U> void
// T and U = stim::simd_bits or stim::simd_bits_range_ref
// Different templates are used because their backing widths
// may differ (e.g. 64 bit vs 256 bit).
copy_where(T from, T to, U pred) {
    from.for_each_word(to, pred, 
            [&] (auto& f, auto& t, auto& p)
            {
                t = p.andnot(f) | (f & p);
            });
}

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

    StateSimulator(const StateSimulator& other)
        :n_qubits(other.n_qubits),
        max_shots(other.max_shots),
        record_table(other.record_table),
        lock_table(other.lock_table),
        record_table_cpy(other.record_table_cpy),
        lock_table_cpy(other.lock_table_cpy),
        rng(other.rng)
    {}

    void    set_seed(uint64_t x) { rng.seed(x); }

    virtual void reset_sim(void) {
        record_table.clear();
        lock_table.clear();
    }
    
    // As operations are expected to be executed on simd_bits, each
    // simulator will have to implement the basic operations themselves.

    // Below: tr = trial. If -1 (default), then that means to execute on all trials.
    // This allows for fine-grained updates. Updates that between all trials and single
    // trials (i.e. a subset of trials) should save the state and rollback after performing
    // a gate, or simply perform the gate for each single trial (a judgement call).
    virtual void    H(std::vector<uint>, int64_t tr=-1) =0;
    virtual void    X(std::vector<uint>, int64_t tr=-1) =0;
    virtual void    Z(std::vector<uint>, int64_t tr=-1) =0;
    virtual void    S(std::vector<uint>, int64_t tr=-1) =0;
    virtual void    CX(std::vector<uint>, int64_t tr=-1) =0;
    virtual void    R(std::vector<uint>, int64_t tr=-1) =0;

    // Other operations that are not guaranteed to have an implementation.
    virtual void    LEAKAGE_ISWAP(std::vector<uint>, int64_t tr=-1) {};

    // As these operations are non-Clifford, we
    // provide a default implementation that does
    // nothing.
    virtual void    T(std::vector<uint>, int64_t tr=-1) {}
    virtual void    RX(fp_t, std::vector<uint>, int64_t tr=-1) {}
    virtual void    RY(fp_t, std::vector<uint>, int64_t tr=-1) {}
    virtual void    RZ(fp_t, std::vector<uint>, int64_t tr=-1) {}

    // Measurement is a special operation. 
    // The first argument is the operands (as with other operations).
    // The second and third arguments (mXwY) is the probability that we
    //  measured a value X but the correct state was Y.
    // The last value is the location in the record table to store the measurement.
    virtual void    M(std::vector<uint>,
                        std::vector<fp_t> m1w0,
                        std::vector<fp_t> m0w1,
                        int record=-1,
                        int64_t tr=-1) =0;
    // 
    // Error operations:
    //      The StateSimulator virtual class will handle rare error injection,
    //      so all the subclass has to do is implement the error phenomenon. 
    //
    // To add new errors, all that must be done is:
    //  (1) Declare the virtual function below (i.e. eXYZ)
    //  (2) Implement the virtual function in the simulator subclasses.
    //      An option is to leave the function blank (as is done with S gates in
    //      FrameSimulator).
    //  + any additional infrastructure you may need.
    //
    typedef void (StateSimulator::*ErrorChannel1Q)(uint, uint64_t);
    typedef void (StateSimulator::*ErrorChannel2Q)(uint, uint, uint64_t);

    template <ErrorChannel1Q CH> void
    error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
        for (size_t i = 0; i < operands.size(); i++) {
            uint j = operands[i];
            stim::RareErrorIterator::for_samples(rates[i], shots, rng,
            [&] (size_t t) {
                if (lock_table[j][t])   return;
                (this->*CH)(j, t);
            });
        }
    }

    template <ErrorChannel2Q CH> void
    error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
        for (size_t i = 0; i < operands.size(); i += 2) {
            uint j1 = operands[i];
            uint j2 = operands[i+1];
            stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
            [&] (size_t t) 
            {
                if (lock_table[j1][t] || lock_table[j2][t]) return;
                (this->*CH)(j1, j2, t);
            });
        }
    }

    virtual void
    error_channel_m(uint64_t rec, fp_t m1w0, fp_t m0w1, stim::simd_bits_range_ref<SIMD_WIDTH> lock) {
        stim::RareErrorIterator::for_samples(m1w0, shots, rng,
                [&] (size_t t)
                {
                    if (lock[t]) return;
                    if (record_table[rec][t] == 0) record_table[rec][t] = 1;
                });
        stim::RareErrorIterator::for_samples(m0w1, shots, rng,
                [&] (size_t t)
                {
                    if (lock[t]) return;
                    if (record_table[rec][t] == 1) record_table[rec][t] = 0;
                });
    }

    virtual void    eDP1(uint, uint64_t) =0;
    virtual void    eX(uint, uint64_t) =0;
    virtual void    eY(uint, uint64_t) =0;
    virtual void    eZ(uint, uint64_t) =0;
    virtual void    eL(uint, uint64_t) =0;

    virtual void    eDP2(uint, uint, uint64_t) =0;
    virtual void    eLI(uint, uint, uint64_t) =0;
    virtual void    eLT(uint, uint, uint64_t) =0;

    void    shift_record_by(uint64_t);
                            // Record shifting is particularly useful for
                            // operating on large programs, where the sliding
                            // window moves over time.

    virtual void    snapshot(void);
                            // Saves the current state of the simulator.
    virtual void rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH>);

    fp_t get_probability_sample_from_rng(void) {
        static std::uniform_real_distribution<> dist(0.0, 1.0);
        return dist(rng);
    }

    stim::simd_bit_table<SIMD_WIDTH>    record_table;
    stim::simd_bit_table<SIMD_WIDTH>    lock_table;
    uint64_t    shots;
protected:
    stim::simd_bit_table<SIMD_WIDTH>    record_table_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    lock_table_cpy;

    std::mt19937_64 rng;

    const uint      n_qubits;
    const uint64_t  max_shots;
};

}   // qontra

#endif  // STATE_SIM_h
