/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#ifndef QONTRA_FRAME_SIM_h
#define QONTRA_FRAME_SIM_h

#include "qontra/sim/base/state_sim.h"

namespace qontra {

class FrameSimulator : public StateSimulator {
public:
    FrameSimulator(uint64_t n, uint64_t max_shots);
    FrameSimulator(const FrameSimulator&);
    FrameSimulator(FrameSimulator&&);

    void reset_sim(void) override;

    void    H(std::vector<uint64_t>, int64_t fr=-1) override;
    void    X(std::vector<uint64_t>, int64_t fr=-1) override;
    void    Z(std::vector<uint64_t>, int64_t fr=-1) override;
    void    S(std::vector<uint64_t>, int64_t fr=-1) override;
    void    CX(std::vector<uint64_t>, int64_t fr=-1) override;
    void    M(std::vector<uint64_t>, std::vector<fp_t>, std::vector<fp_t>, int record=-1, int64_t fr=-1) override;
    void    R(std::vector<uint64_t>, int64_t fr=-1) override;

    void    LEAKAGE_ISWAP(std::vector<uint64_t>, int64_t fr=-1) override;

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
private:
    stim::simd_bit_table<SIMD_WIDTH>    x_table;
    stim::simd_bit_table<SIMD_WIDTH>    z_table;
    stim::simd_bit_table<SIMD_WIDTH>    leak_table;

    stim::simd_bit_table<SIMD_WIDTH>    x_table_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    z_table_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    leak_table_cpy;

    friend class MemorySimulator;
};

}   // qontra

#include "frame_sim.inl"

#endif  // QONTRA_FRAME_SIM_h
