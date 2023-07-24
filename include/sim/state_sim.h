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
#include <tuple>
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

    virtual void reset_sim(void) {
        record_table.clear();
        lock_table.clear();

        no_error_nlog_probability = 0.0;
    }
    
    // As operations are expected to be executed on simd_bits, each
    // simulator will have to implement the basic operations themselves.

    virtual void    H(std::vector<uint>) =0;
    virtual void    X(std::vector<uint>) =0;
    virtual void    Z(std::vector<uint>) =0;
    virtual void    S(std::vector<uint>) =0;
    virtual void    CX(std::vector<uint>) =0;
    virtual void    M(std::vector<uint>, int record=-1) =0;
    virtual void    R(std::vector<uint>) =0;

    enum class ErrorLabel { I, X, Y, Z, L };
    // 
    // Error operations:
    //      The StateSimulator virtual class will handle rare error injection,
    //      so all the subclass has to do is implement the error phenomenon. 
    //
    //      Each error function must return an ErrorLabel (or a pair of ErrorLabels)
    //      to indicate how each qubit was affected. That is,
    //          I = qubit had no error
    //          X = qubit had an X error
    //          Y = qubit had a Y error
    //          Z = qubit had a Z error
    //          L = qubit had a leakage error
    //      Y covers the case where X and Z both occur. L is mutually exclusive
    //      compared to the other errors (as a leaked qubit does not care about
    //      other errors). The ErrorLabels are used to build an error log that
    //      can be used to record exactly what errors have occurred.
    //
    //      The types of errors possible can be extended, but we limit the labels
    //      to Pauli+leakage as these are standard. One possible extension is
    //      adding errors such as rt(X), rt(Y), rt(Z) and their inverses.
    //
    // To add new errors, all that must be done is:
    //  (1) Declare the virtual function below (i.e. eXYZ)
    //  (2) Implement the virtual function in the simulator subclasses.
    //      An option is to leave the function blank (as is done with S gates in
    //      FrameSimulator).
    //  + any additional infrastructure you may need.
    //
    // K-qubit operations can be implemented by adding a new error channel function
    // along with a new typedef (i.e. label_kq_t). We do not natively provide this
    // as 2-qubit gates are standard on most platforms.
    //      
    typedef ErrorLabel                          label_1q_t;
    typedef std::pair<ErrorLabel, ErrorLabel>   label_2q_t;
    typedef label_1q_t (StateSimulator::*ErrorChannel1Q)(uint, uint64_t);
    typedef label_2q_t (StateSimulator::*ErrorChannel2Q)(uint, uint, uint64_t);

    template <ErrorChannel1Q CH> void 
    error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
        for (uint i = 0; i < operands.size(); i++) {
            uint j = operands[i];
            stim::RareErrorIterator::for_samples(rates[i], shots, rng,
                    [&] (size_t t) {
                        if (lock_table[j][t])   return;
                        label_1q_t label = (this->*CH)(j, t);
                        error_log.push_back({
                                (int)j, 
                                -1, 
                                t, 
                                rates[i], 
                                label,
                                ErrorLabel::I
                        });
                    });
            no_error_nlog_probability += -log(1 - rates[i]);
        }
    }

    template <ErrorChannel2Q CH> void
    error_channel(std::vector<uint> operands, std::vector<fp_t> rates) {
        for (uint i = 0; i < operands.size(); i += 2) {
            uint j1 = operands[i];
            uint j2 = operands[i+1];
            stim::RareErrorIterator::for_samples(rates[i>>1], shots, rng,
                    [&] (size_t t) 
                    {
                        if (lock_table[j1][t] || lock_table[j2][t]) return;
                        label_2q_t labels = (this->*CH)(j1, j2, t);
                        error_log.push_back({
                            (int)j1, 
                            (int)j2, 
                            t,
                            rates[i>>1],
                            labels.first,
                            labels.second
                        });
                    });
            no_error_nlog_probability += -log(1 - rates[i>>1]);
        }
    }

    virtual label_1q_t  eDP1(uint, uint64_t) =0;
    virtual label_1q_t  eX(uint, uint64_t) =0;
    virtual label_1q_t  eY(uint, uint64_t) =0;
    virtual label_1q_t  eZ(uint, uint64_t) =0;
    virtual label_1q_t  eL(uint, uint64_t) =0;

    virtual label_2q_t  eDP2(uint, uint, uint64_t) =0;
    virtual label_2q_t  eLI(uint, uint, uint64_t) =0;
    virtual label_2q_t  eLT(uint, uint, uint64_t) =0;

    void    shift_record_by(uint64_t);
                            // Record shifting is particularly useful for
                            // operating on large programs, where the sliding
                            // window moves over time.

    virtual void    snapshot(void);
                            // Saves the current state of the simulator.
    virtual void    rollback_where(stim::simd_bits_range_ref);
                            // Rolls back the state to the snapshot

    typedef struct {
        int         q1 = -1;    // A qubit that is <0 is invalid.
        int         q2 = -1;
        uint64_t    trial;

        fp_t    error_probability;

        ErrorLabel  error_on_q1 = ErrorLabel::I;
        ErrorLabel  error_on_q2 = ErrorLabel::I;
    } error_data_t;

    std::vector<error_data_t>   get_error_log(void) { return error_log; }
    void                        clear_error_log(void) { error_log.clear(); }

    fp_t get_negative_log_probability_of_no_error(void) {
        return no_error_nlog_probability;
    }

    stim::simd_bit_table    record_table;
    stim::simd_bit_table    lock_table;
    uint64_t                shots;
protected:
    stim::simd_bit_table    record_table_cpy;
    stim::simd_bit_table    lock_table_cpy;

    std::mt19937_64 rng;

    const uint      n_qubits;
    const uint64_t  max_shots;
private:
    std::vector<error_data_t>   error_log;

    fp_t    no_error_nlog_probability;  // We track the probability that no
                                        // error will occur. We use -log
                                        // probability to avoid floating point
                                        // error.
};

void    copy_where(stim::simd_bits_range_ref from,
                    stim::simd_bits_range_ref to,
                    stim::simd_bits_range_ref pred);

}   // qontra

#endif  // STATE_SIM_h
