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
            apply_1q_gate(h_base(x));
        }
    }

    void X(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(x_base(x));
        }
    }

    void Z(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(z_base(x));
        }
    }

    void S(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(s_base(x));
        }
    }

    void T(std::vector<uint> operands) override {
        for (uint x : operands) {
            apply_1q_gate(t_base(x));
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
            apply_1q_gate(rx_base(angle, x), x);
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

    void    M(std::vector<uint>, fp_t, fp_t, int record=-1) override;
    void    R(std::vector<uint>) override;

    StateSimulator::label_1q_t  eDP1(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eX(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eY(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eZ(uint, uint64_t=0) override;
    StateSimulator::label_1q_t  eL(uint, uint64_t=0) override;

    StateSimulator::label_2q_t  eDP2(uint, uint, uint64_t=0) override;
    StateSimulator::label_2q_t  eLI(uint, uint, uint64_t=0) override;
    StateSimulator::label_2q_t  eLT(uint, uint, uint64_t=0) override;

    void    snapshot(void) override;
    void    rollback_where(stim::simd_bits_range_ref) override;
    
    struct {
        // Simulation parameters.
        fp_t    leakage_transport_angle = 0.2*M_PI; // RY angle of leakage transport 
                                                    //  (0 = no transport, PI = always transport).
        fp_t    faulty_cx_angle = 0.65*M_PI;        // RX error caused by leaked state.

        // If any of these are true, then the corresponding value
        // is taken from the above variables. Otherwise, the angle
        // is randomized.
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

    cytnx::UniTensor    build_1q_gate(const cytnx::Tensor&, uint);
    cytnx::UniTensor    build_2q_gate(const cytnx::Tensor&, uint, uint);
    void                apply_1q_gate(const cytnx::UniTensor, uint);
    void                apply_2q_gate(const cytnx::UniTensor, uint, uint);

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
