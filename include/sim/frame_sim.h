/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#ifndef FRAME_SIM_h
#define FRAME_SIM_h

#include "sim/state_sim.h"

#include <vector>

namespace qontra {

class FrameSimulator : public StateSimulator {
public:
    FrameSimulator(uint n, uint64_t max_shots);
    FrameSimulator(const FrameSimulator&)
    FrameSimulator(FrameSimulator&&);

    void reset_sim(void) override;

    void    H(std::vector<uint>, int64_t fr=-1) override;
    void    X(std::vector<uint>, int64_t fr=-1) override;
    void    Z(std::vector<uint>, int64_t fr=-1) override;
    void    S(std::vector<uint>, int64_t fr=-1) override;
    void    CX(std::vector<uint>, int64_t fr=-1) override;
    void    M(std::vector<uint>, std::vector<fp_t>, std::vector<fp_t>, int record=-1, int64_t fr=-1) override;
    void    R(std::vector<uint>, int64_t fr=-1) override;

    void    LEAKAGE_ISWAP(std::vector<uint>, int64_t fr=-1) override;

    void    eDP1(uint, uint64_t) override;
    void    eX(uint, uint64_t) override;
    void    eY(uint, uint64_t) override;
    void    eZ(uint, uint64_t) override;
    void    eL(uint, uint64_t) override;

    void    eDP2(uint, uint, uint64_t) override;
    void    eLI(uint, uint, uint64_t) override;
    void    eLT(uint, uint, uint64_t) override;

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

#endif  // FRAME_SIM_h
