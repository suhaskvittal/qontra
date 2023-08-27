/*
 *  author: Suhas Vittal
 *  date:   21 August 2023
 *
 *  This is a tensor network simulator.
 *  Unlike the other simulators,
 *  this is only compatible with
 *  G_SHOTS_PER_BATCH = 1.
 * */

#ifndef UNIVERSAL_SIM_h
#define UNIVERSAL_SIM_h

#include "experiments.h"
#include "sim/state_sim.h"

#include <UniTensor.hpp>

#include <complex>
#include <random>
#include <string>
#include <vector>

#include <assert.h>
#include <math.h>

namespace qontra {

template <int W=2>
class UniversalSimulator : StateSimulator {
public:
    UniversalSimulator(uint n);

    void    reset_sim(void) override;

    void H(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(h_base(x), x);
        }
    }

    void X(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(x_base(x), x);
        }
    }

    void Z(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(z_base(x), x);
        }
    }

    void S(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(s_base(x), x);
        }
    }

    void T(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(t_base(x), x);
        }
    }

    void CX(std::vector<uint> operands) override {
        for (uint i = 0; i < operands.size(); i += 2) {
            uint x = operands[i], y = operands[i+1];
            apply_2q_gate(cx_base(x, y), x, y);
        }
    }
    void RX(fp_t angle, std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(override(angle, x), x);
        }
    }

    void RY(fp_t angle, std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(ry_base(angle, x), x);
        }
    }

    void RZ(fp_t angle, std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(rz_base(angle, x), x);
        }
    }

    void R(std::vector<uint> operands) override {
        for (uint x : operands) {
            auto tmp = test_1q_gate(reset_base(x), x);
            tmp.normalize_();
            update_1q_state(tmp, x);
        }
    }

    void    M(std::vector<uint>, fp_t, fp_t, int record=-1) override;

    StateSimulator::label_1q_t  eDP1(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eX(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eY(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eZ(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eL(uint, uint64_t=0) override;

    // As we are performing a tensor network simulation, we
    // might as well implement an amplitude-damping and phase
    // damping channel. However, a good implementation of this
    // channel, as in Google's Sycamore experiments, requires
    // a density matrix. However, density matrices are not good
    // for non-determinism, which we seek to emulate. So, we create
    // three sub-channels:
    //
    // Amplitude Damping (AD):
    //  (1) relaxes |x+1> --> |x> via a random X rotation
    // Phase Damping (PD):
    //  (3) an RZ rotation in the |x>-|x+1> basis
    // Heating (HEAT):
    //  (4) excites |x> --> |x+1> via an X rotation
    //
    // The b1, b2 template parameters perform the channels in
    // the |b1>-|b2> basis.
    template <int b1, int b2>   StateSimulator::label_1q_t  eAD(uint, uint64_t=0);
    template <int b1, int b2>   StateSimulator::label_1q_t  ePD(uint, uint64_t=0);
    template <int b1, int b2>   StateSimulator::label_1q_t  eHEAT(uint, uint64_t=0);

    StateSimulator::label_2q_t  eDP2(uint, uint, uint64_t=0) override;
    StateSimulator::label_2q_t  eLI(uint, uint, uint64_t=0) override;
    StateSimulator::label_2q_t  eLT(uint, uint, uint64_t=0) override;

    void    snapshot(void) override;
    void    rollback_where(stim::simd_bits_range_ref) override;
    
    struct {
        // Simulation parameters.
        fp_t    amplitude_damping_angle = 0.85*M_PI;    // RX angle of amplitude damping
        fp_t    phase_damping_angle = 0.85*M_PI;        // RZ angle of phase damping
        fp_t    heating_angle = 0.6*M_PI;               // RX angle of heating excitation
        fp_t    leakage_transport_angle = 0.2*M_PI;     // RY angle of leakage transport 
                                                        //  (0 = no transport, PI = always transport).
        fp_t    faulty_cx_angle = 0.65*M_PI;            // RX error caused by leaked state.

        // If any of these are true, then the corresponding value
        // is taken from the above variables. Otherwise, the angle
        // is randomized.
        bool    fix_amplitude_damping_angle = true;
        bool    fix_phase_damping_angle = true;
        bool    fix_heating_angle = true;
        bool    fix_leakage_transport_angle = true;
        bool    fix_faulty_cx_angle = false;
    } params;
private:
    cytnx::UniTensor   qubit_base(int);

    cytnx::UniTensor    h_base(uint);
    cytnx::UniTensor    x_base(uint);
    cytnx::UniTensor    y_base(uint);
    cytnx::UniTensor    z_base(uint);
    cytnx::UniTensor    s_base(uint);
    cytnx::UniTensor    t_base(uint);
    cytnx::UniTensor    cx_base(uint, uint);
    cytnx::UniTensor    rx_base(fp_t, uint);
    cytnx::UniTensor    ry_base(fp_t, uint);
    cytnx::UniTensor    rz_base(fp_t, uint);

    cytnx::UniTensor    m_base(uint, uint b);

    cytnx::UniTensor    reset_base(uint);

    cytnx::UniTensor    build_1q_gate(const cytnx::Tensor&, uint);
    cytnx::UniTensor    build_2q_gate(const cytnx::Tensor&, uint, uint);
    cytnx::UniTensor    test_1q_gate(const cytnx::UniTensor, uint);
    cytnx::UniTensor    test_2q_gate(const cytnx::UniTensor, uint, uint);

    void    apply_1q_gate(const cytnx::UniTensor, uint);
    void    update_1q_state(const cytnx::UniTensor, uint);
    void    apply_2q_gate(const cytnx::UniTensor, uint, uint);
    void    update_2q_state(const cytnx::UniTensor, uint, uint);

    bool    is_bond_length_critical(void);
    void    compress(void);

    std::map<uint, std::string>                 qubit_to_label;
    std::set<cytnx::UniTensor*>                 tensor_pool;
    std::map<std::string, cytnx::UniTensor*>    label_to_tensor;
    std::map<cytnx::UniTensor*, std::vector<std::string>>   tensor_to_labels;
    std::map<uint, std::string> awaiting_labels;    // New labels after operation.

    int32_t label_ctr;

    std::uniform_real_distribution  rpdist;
};

// build_leakage_operator and build_rotation_operator_for_b2_subspace
// are used to construct the leakage operator for control leakage in this
// paper: 
//  "Efficient Simulation of Leakage Errors in
//  Quantum Error Correcting Codes Using Tensor Network Methods"
//      by Hidetaka Manabe, Yasunari Suzuki, Andrew S. Darmawan
//  See Section D.

cytnx::Tensor
build_leakage_operator(void);

template <int b> cyntx::Tensor
build_rotation_operator_for_b2_subspace(fp_t, fp_t);


}   // qontra

#endif  // UNIVERSAL_SIM_h
