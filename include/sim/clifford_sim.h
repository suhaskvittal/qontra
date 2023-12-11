/* author: Suhas Vittal
 *  date:   18 June 2023
 *
 *  Fast Clifford simulator.
 * */

#ifndef CLIFFORD_SIM_h
#define CLIFFORD_SIM_h

#include "sim/state_sim.h"

#include <vector>

#include <stim.h>

namespace qontra {

#define __CSIM_X_WIDTH(n)   2*n*(n+1)+1
#define __CSIM_Z_WIDTH(n)   2*n*(n+1)
#define __CSIM_R_WIDTH(n)   2*n+4

class CliffordSimulator : public StateSimulator {
public:
    CliffordSimulator(uint n, uint64_t max_shots)
        :StateSimulator(n, max_shots),
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
        leak_table_cpy(n, max_shots)
    {
        init_tables();
    }

    CliffordSimulator(const CliffordSimulator& other)
        :StateSimulator(other),
        x_width(other.x_width),
        z_width(other.z_width),
        r_width(other.r_width),
        leak_width(other.leak_width),
        x_table(other.x_table),
        z_table(other.z_table),
        r_table(other.r_table),
        leak_table(other.leak_table),
        x_table_cpy(other.x_table_cpy),
        z_table_cpy(other.z_table_cpy),
        r_table_cpy(other.r_table_cpy),
        leak_table_cpy(other.leak_table_cpy)
    {}

    void reset_sim(void) override {
        StateSimulator::reset_sim();

        x_table.clear();
        z_table.clear();
        r_table.clear();
        leak_table.clear();
        init_tables();
    }

    void    H(std::vector<uint>) override;
    void    X(std::vector<uint>) override;
    void    Z(std::vector<uint>) override;
    void    S(std::vector<uint>) override;
    void    CX(std::vector<uint>) override;
    void    M(std::vector<uint>, std::vector<fp_t>, std::vector<fp_t>, int record=-1) override;
    void    R(std::vector<uint>) override;

    StateSimulator::label_1q_t  eDP1(uint, uint64_t) override;
    StateSimulator::label_1q_t  eX(uint, uint64_t) override;
    StateSimulator::label_1q_t  eY(uint, uint64_t) override;
    StateSimulator::label_1q_t  eZ(uint, uint64_t) override;
    StateSimulator::label_1q_t  eL(uint, uint64_t) override;

    StateSimulator::label_2q_t  eDP2(uint, uint, uint64_t) override;
    StateSimulator::label_2q_t  eLI(uint, uint, uint64_t) override;
    StateSimulator::label_2q_t  eLT(uint, uint, uint64_t) override;

    void    snapshot(void) override;
    void    rollback_where(stim::simd_bits_range_ref) override;
private:
    void    init_tables(void);
    void    rowsum(uint h, uint i, uint64_t t); // Only executes rowsum on one trial.
    void    browsum(uint h, uint i, bool use_pred, stim::simd_bits_range_ref pred);
                // Performs rowsum on all trials at the same time, and 
                // conditions the update based on the predicate.

    stim::simd_bit_table    x_table;    // Has an extra row for scratch space.
    stim::simd_bit_table    z_table;
    stim::simd_bit_table    r_table;    // Has four extra rows for executing rowsum.
    stim::simd_bit_table    leak_table; // Only has n rows.

    stim::simd_bit_table    x_table_cpy;
    stim::simd_bit_table    z_table_cpy;
    stim::simd_bit_table    r_table_cpy;
    stim::simd_bit_table    leak_table_cpy;

    const uint x_width;
    const uint z_width;
    const uint r_width;
    const uint leak_width;
};

}   // qontra

#endif  // CLIFFORD_SIM_h
