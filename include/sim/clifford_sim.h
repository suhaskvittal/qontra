/*
 *  author: Suhas Vittal
 *  date:   18 June 2023
 *
 *  Fast Clifford simulator.
 * */

#ifndef CLIFFORD_SIM_h
#define CLIFFORD_SIM_h

#include "defs.h"

#include <random>
#include <vector>

#include <stim.h>

namespace qontra {

#define __CSIM_X_WIDTH(n)   2*n*(n+1)+1
#define __CSIM_Z_WIDTH(n)   2*n*(n+1)
#define __CSIM_R_WIDTH(n)   2*n+4

namespace cliffsim {

extern uint64_t G_RECORD_SPACE_SIZE;

}   //  cliffsim

class CliffordSimulator {
public:
    CliffordSimulator(uint n, uint64_t max_shots)
        :n_qubits(n),
        max_shots(max_shots),
        x_width(__CSIM_X_WIDTH(n)),
        z_width(__CSIM_Z_WIDTH(n)),
        r_width(__CSIM_R_WIDTH(n)),
        leak_width(n),
        x_table(__CSIM_X_WIDTH(n), max_shots),
        z_table(__CSIM_Z_WIDTH(n), max_shots),
        r_table(__CSIM_R_WIDTH(n), max_shots),
        leak_table(n, max_shots),
        x_table_cpy(__CSIM_X_WIDTH(n), max_shots),
        z_table_cpy(__CSIM_Z_WIDTH(n), max_shots),
        r_table_cpy(__CSIM_R_WIDTH(n), max_shots),
        leak_table_cpy(n, max_shots),
        record_table(cliffsim::G_RECORD_SPACE_SIZE, max_shots),
        record_offset(0),
        rng(0)
    {
        init_tables();
    }

    void reset_sim(void) {
        x_table.clear();
        z_table.clear();
        r_table.clear();
        leak_table.clear();
        record_table.clear();

        record_offset = 0;

        init_tables();
    }

    void set_seed(uint x)   { rng.seed(x); }

    void    H(std::vector<uint>);
    void    X(std::vector<uint>);
    void    Z(std::vector<uint>);
    void    S(std::vector<uint>);
    void    CX(std::vector<uint>);
    void    M(std::vector<uint>, bool record=true);
    void    R(std::vector<uint>);

    void    eDPO(std::vector<uint>, std::vector<fp_t>);
    void    eDPH(std::vector<uint>, std::vector<fp_t>);

    void    eDP1(std::vector<uint>, std::vector<fp_t>);
    void    eDP2(std::vector<uint>, std::vector<fp_t>);
    void    eX(std::vector<uint>, std::vector<fp_t>);
    void    eLI(std::vector<uint>, std::vector<fp_t>);
    void    eLT(std::vector<uint>, std::vector<fp_t>);

    void    reduce_record_by(uint64_t);
    void    shift_record_by(uint64_t);

    void    snapshot(void); // Copies the state tables into new data structures
                            // for use later.
    void    rollback_at_trial(uint64_t);    // Rolls back the state to the snapshot

    uint64_t    get_record_size(void) { return record_offset; }

    stim::simd_bit_table    record_table;
    uint64_t                shots;
private:
    void    init_tables(void);
    void    rowsum(uint h, uint i, bool use_pred, stim::simd_bits_range_ref pred);
                // Conditions the rowsum based on the predicate.

    stim::simd_bit_table    x_table;    // Has an extra row for scratch space.
    stim::simd_bit_table    z_table;
    stim::simd_bit_table    r_table;    // Has four extra rows for executing rowsum.
    stim::simd_bit_table    leak_table; // Only has n rows.

    stim::simd_bit_table    x_table_cpy;
    stim::simd_bit_table    z_table_cpy;
    stim::simd_bit_table    r_table_cpy;
    stim::simd_bit_table    leak_table_cpy;

    const uint      n_qubits;
    const uint64_t  max_shots;

    const uint x_width;
    const uint z_width;
    const uint r_width;
    const uint leak_width;

    uint64_t record_offset;

    std::mt19937_64 rng;
};

}   // qontra

#endif  // CLIFFORD_SIM_h
