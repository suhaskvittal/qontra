/* author: Suhas Vittal
 *  date:   18 June 2023
 *
 *  Fast Clifford simulator.
 * */

#ifndef QONTRA_CLIFFORD_SIM_h
#define QONTRA_CLIFFORD_SIM_h

#include "qontra/sim/base/state_sim.h"

namespace qontra {

class CliffordSimulator : public StateSimulator {
public:
    CliffordSimulator(uint64_t n, uint64_t max_shots);
    CliffordSimulator(const CliffordSimulator& other);

    void    reset_sim(void) override;

    void    H(std::vector<uint64_t>, int64_t tr=-1) override;
    void    X(std::vector<uint64_t>, int64_t tr=-1) override;
    void    Z(std::vector<uint64_t>, int64_t tr=-1) override;
    void    S(std::vector<uint64_t>, int64_t tr=-1) override;
    void    CX(std::vector<uint64_t>, int64_t tr=-1) override;
    void    M(std::vector<uint64_t>, std::vector<fp_t>, std::vector<fp_t>, int record=-1, int64_t tr=-1) override;
    void    R(std::vector<uint64_t>, int64_t tr=-1) override;

    void    eDP1(uint64_t, uint64_t) override;
    void    eX(uint64_t, uint64_t) override;
    void    eY(uint64_t, uint64_t) override;
    void    eZ(uint64_t, uint64_t) override;
    void    eL(uint64_t, uint64_t) override;

    void    eDP2(uint64_t, uint64_t, uint64_t) override;
    void    eLI(uint64_t, uint64_t, uint64_t) override;
    void    eLT(uint64_t, uint64_t, uint64_t) override;

    void    snapshot(void) override;
    void    rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH>) override;

    stim::simd_bit_table<SIMD_WIDTH>    x_table;    // Has an extra row for scratch space.
    stim::simd_bit_table<SIMD_WIDTH>    z_table;
    stim::simd_bit_table<SIMD_WIDTH>    r_table;    // Has four extra rows for executing rowsum.
private:
    uint64_t  get_index(uint64_t i, uint64_t j);

    void    clear_scratch_space(void);
    void    clear_row(uint64_t);
    void    clear_row_for_trial(uint64_t, int64_t);
    void    clear_row_where(uint64_t, stim::simd_bits_range_ref<SIMD_WIDTH>);

    void    swap_rows_for_trial(uint64_t, uint64_t, int64_t);
    void    swap_rows_where(uint64_t, uint64_t, stim::simd_bits_range_ref<SIMD_WIDTH>);

    void    batch_measure_qubit(uint64_t);
    void    measure_qubit_in_trial(uint64_t, int64_t);
        
    void    init_tables(void);

    void    rowsum(uint64_t h, uint64_t i, uint64_t t);
                // Only executes rowsum on one trial.
    void    browsum(uint64_t h, uint64_t i, stim::simd_bits_range_ref<SIMD_WIDTH> pred);
                // Performs rowsum on all trials at the same time, and 
                // conditions the update based on the predicate.

    stim::simd_bit_table<SIMD_WIDTH>    x_table_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    z_table_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    r_table_cpy;

    const uint64_t x_width;
    const uint64_t z_width;
    const uint64_t r_width;
};

}   // qontra

#include "clifford_sim.inl"

#endif  // QONTRA_CLIFFORD_SIM_h
