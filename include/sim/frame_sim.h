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
    FrameSimulator(uint n, uint64_t max_shots)
        :StateSimulator(n, max_shots),
        x_table(n, max_shots),
        z_table(n, max_shots),
        leak_table(n, max_shots),
        x_table_cpy(n, max_shots),
        z_table_cpy(n, max_shots),
        leak_table_cpy(n, max_shots)
    {}

    FrameSimulator(const FrameSimulator& other)
        :StateSimulator(other),
        x_table(other.x_table),
        z_table(other.z_table),
        leak_table(other.leak_table),
        x_table_cpy(other.x_table_cpy),
        z_table_cpy(other.z_table_cpy),
        leak_table_cpy(other.leak_table_cpy)
    {}

    void reset_sim(void) override {
        StateSimulator::reset_sim();

        x_table.clear();
        z_table.clear();
        leak_table.clear();
    }

    void    H(std::vector<uint>, int64_t fr=-1) override;
    void    X(std::vector<uint>, int64_t fr=-1) override;
    void    Z(std::vector<uint>, int64_t fr=-1) override;
    void    S(std::vector<uint>, int64_t fr=-1) override;
    void    CX(std::vector<uint>, int64_t fr=-1) override;
    void    M(std::vector<uint>, std::vector<fp_t>, std::vector<fp_t>, int record=-1, int64_t fr=-1) override;
    void    R(std::vector<uint>, int64_t fr=-1) override;

    void    LEAKAGE_ISWAP(std::vector<uint>, int64_t fr=-1) override;

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
    stim::simd_bit_table    x_table;
    stim::simd_bit_table    z_table;
    stim::simd_bit_table    leak_table;

    stim::simd_bit_table    x_table_cpy;
    stim::simd_bit_table    z_table_cpy;
    stim::simd_bit_table    leak_table_cpy;

    friend class MemorySimulator;
};

}   // qontra

#endif  // FRAME_SIM_h
