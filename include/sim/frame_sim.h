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

    void reset_sim(void) override {
        StateSimulator::reset_sim();

        x_table.clear();
        z_table.clear();
        leak_table.clear();
    }

    void    H(std::vector<uint>) override;
    void    X(std::vector<uint>) override;
    void    Z(std::vector<uint>) override;
    void    S(std::vector<uint>) override;
    void    CX(std::vector<uint>) override;
    void    M(std::vector<uint>, bool record=true) override;
    void    R(std::vector<uint>) override;

    void    eDPO(std::vector<uint>, std::vector<fp_t>) override;
    void    eDPH(std::vector<uint>, std::vector<fp_t>) override;

    void    eDP1(std::vector<uint>, std::vector<fp_t>) override;
    void    eDP2(std::vector<uint>, std::vector<fp_t>) override;
    void    eX(std::vector<uint>, std::vector<fp_t>) override;
    void    eLI(std::vector<uint>, std::vector<fp_t>) override;
    void    eLT(std::vector<uint>, std::vector<fp_t>) override;

    void    snapshot(void) override;
    void    rollback_at_trial(uint64_t) override;
private:
    stim::simd_bit_table    x_table;
    stim::simd_bit_table    z_table;
    stim::simd_bit_table    leak_table;

    stim::simd_bit_table    x_table_cpy;
    stim::simd_bit_table    z_table_cpy;
    stim::simd_bit_table    leak_table_cpy;
};

}   // qontra

#endif  // FRAME_SIM_h
