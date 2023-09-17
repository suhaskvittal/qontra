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

#include <cytnx.hpp>
#include <Scalar.hpp>
#include <Tensor.hpp>
#include <UniTensor.hpp>

#include <complex>
#include <random>
#include <string>
#include <vector>

#include <assert.h>
#include <math.h>

namespace qontra {

typedef std::complex<fp_t>  c64;
const c64 M_J(0, 1);

template <int W=2>
class UniversalSimulator : public StateSimulator {
public:
    UniversalSimulator(uint n)
        :StateSimulator(n, 1),
        qubit_to_label(),
        tensor_pool(),
        label_to_tensor(),
        tensor_to_labels(),
        awaiting_labels(),
        label_ctr(n+1),
        rpdist(0.0, 1.0)
    {
        assert(experiments::G_SHOTS_PER_BATCH == 1);
        reset_sim();
    }

    ~UniversalSimulator() {
        for (cytnx::UniTensor* x : tensor_pool)   delete x;
    }

    void reset_sim(void) override {
        // Delete all tensors in the pool.
        for (cytnx::UniTensor* x : tensor_pool)   delete x;
        tensor_pool.clear();
        label_to_tensor.clear();
        tensor_to_labels.clear();
        // Re-instantiate the starting tensors.
        for (int i = 0; i < n_qubits; i++) {
            std::string label = std::to_string(i+1);
            cytnx::UniTensor* init = new cytnx::UniTensor(qubit_base(i));

            qubit_to_label[i] = label;
            label_to_tensor[label] = init;
            tensor_to_labels[init].push_back(label);
            tensor_pool.insert(init);
        }
    }

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

    void R(std::vector<uint> operands) override {
        for (uint x : operands) {
            auto tmp = test_1q_gate(reset_base(x), x);
            tmp.normalize_();
            update_1q_state(tmp, x);
        }
    }

    void M(std::vector<uint> operands, std::vector<fp_t> m1w0, std::vector<fp_t> m0w1, int record=-1) override {
        uint k = 0;
        for (uint x : operands) {
            fp_t r = rpdist(rng);
            // Get the probability of measuring each state.
            for (uint i = 0; i < W; i++) {
                auto tmp = test_1q_gate(m_base(x, i), x);
                fp_t pr = (fp_t) cytnx::Scalar(tmp.Norm().at({0})).abs();
                pr *= pr;
                if (r < pr) {
                    // Collapse to this state.
                    tmp.normalize_();
                    update_1q_state(tmp, x);
                    if (record >= 0) {
                        uint m = i;
                        // Randomly resolve the measurement if we collapsed to |2>, |3>, etc.
                        if (m > 1) {
                            fp_t x = rpdist(rng);
                            m = x < 0.5 ? 0 : 1;
                        }
                        // Apply errors.
                        fp_t e = rpdist(rng);
                        if (m == 0 && e < m1w0[k])         m = 1;
                        else if (m == 1 && e < m0w1[k])    m = 0;

                        record_table[record++][0] = m;
                    }
                    k++;
                    break;
                } else {
                    r -= pr;
                }
            }
        }
    }

    StateSimulator::label_1q_t eDP1(uint x, uint64_t t=0) override {
        const StateSimulator::label_1q_t IXYZ[4] = { 
            StateSimulator::label_1q_t::I,
            StateSimulator::label_1q_t::X,
            StateSimulator::label_1q_t::Y,
            StateSimulator::label_1q_t::Z
        };

        auto p = rng() & 3;
        if (p & 1)  apply_1q_gate(x_base(x), x);
        if (p & 2)  apply_1q_gate(z_base(x), x);
        return IXYZ[p];
    }

    StateSimulator::label_1q_t eX(uint x, uint64_t t=0) override {
        apply_1q_gate(x_base(x), x);
        return StateSimulator::label_1q_t::X;
    }

    StateSimulator::label_1q_t eY(uint x, uint64_t t=0) override {
        apply_1q_gate(y_base(x), x);
        return StateSimulator::label_1q_t::Y;
    }

    StateSimulator::label_1q_t eZ(uint x, uint64_t t=0) override {
        apply_1q_gate(z_base(x), x);
        return StateSimulator::label_1q_t::Z;
    }

    StateSimulator::label_1q_t eL(uint x, uint64_t t=0) override {
        // Leave this for explicitly specialized implementations.
        return StateSimulator::label_1q_t::I;
    }

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
    template <int b1, int b2> StateSimulator::label_1q_t
    eAD(uint x, uint64_t t=0) {
        fp_t angle;
        if (params.fix_amplitude_damping_angle) {
            angle = params.amplitude_damping_angle;
        } else {
            angle = rpdist(rng) * 2 * M_PI;
        }

        cytnx::Tensor rot({W, W}, cytnx::Type.ComplexDouble);
        for (int i = 0; i < W; i++) {
            if (i != b1 && i != b2) rot.at({i, i}) = 1;
        }
        rot.at({b1, b1}) = 1.0;
        rot.at({b1, b2}) = c64(0.0, -sin(angle*0.5));
        rot.at({b2, b2}) = cos(angle*0.5);

        apply_1q_gate(build_1q_gate(rot, x), x);
        return b1 <= 1 && b2 <= 1 ? StateSimulator::label_1q_t::X : StateSimulator::label_1q_t::L;
    }

    template <int b1, int b2> StateSimulator::label_1q_t
    ePD(uint x, uint64_t t=0) {
        fp_t angle;
        if (params.fix_phase_damping_angle) {
            angle = params.phase_damping_angle;
        } else {
            angle = rpdist(rng) * 2 * M_PI;
        }

        cytnx::Tensor rot({W, W}, cytnx::Type.ComplexDouble);
        for (int i = 0; i < W; i++) {
            if (i != b1 && i != b2) rot.at({i, i}) = 1;
        }
        rot.at({b1, b1}) = c64(cos(angle*0.5), -sin(angle*0.5));
        rot.at({b2, b2}) = c64(cos(angle*0.5), sin(angle*0.5));
        
        apply_1q_gate(build_1q_gate(rot, x), x);
        return b1 <= 1 && b2 <= 1 ? StateSimulator::label_1q_t::Z : StateSimulator::label_1q_t::I;
    }

    template <int b1, int b2> StateSimulator::label_1q_t 
    eHEAT(uint x, uint64_t t=0) {
        fp_t angle;
        if (params.fix_heating_angle) {
            angle = params.heating_angle;
        } else {
            angle = rpdist(rng) * 2 * M_PI;
        }

        cytnx::Tensor rot({W, W}, cytnx::Type.ComplexDouble);
        for (int i = 0; i < W; i++) {
            if (i != b1 && i != b2) rot.at({i, i}) = 1;
        }
        rot.at({b1, b1}) = cos(angle*0.5);
        rot.at({b2, b1}) = c64(0.0, -sin(angle*0.5));
        rot.at({b2, b2}) = 1.0;

        apply_1q_gate(build_1q_gate(rot, x), x);
        return b1 <= 1 && b2 <= 1 ? StateSimulator::label_1q_t::X : StateSimulator::label_1q_t::L;
    }

    StateSimulator::label_2q_t eDP2(uint x, uint y, uint64_t t=0) override {
        const StateSimulator::label_1q_t IXYZ[4] = { 
            StateSimulator::label_1q_t::I,
            StateSimulator::label_1q_t::X,
            StateSimulator::label_1q_t::Y,
            StateSimulator::label_1q_t::Z
        };

        auto p = rng() & 15;
        if (p & 1)  apply_1q_gate(x_base(x), x);
        if (p & 2)  apply_1q_gate(z_base(x), x);
        if (p & 4)  apply_1q_gate(x_base(y), y);
        if (p & 8)  apply_1q_gate(z_base(y), y);
        return std::make_pair(IXYZ[p & 3], IXYZ[(p>>2) & 3]);
    }

    StateSimulator::label_2q_t eLI(uint, uint, uint64_t=0) override {
        const StateSimulator::label_1q_t L_I = StateSimulator::label_1q_t::I;
        return std::make_pair(L_I, L_I);
    }

    StateSimulator::label_2q_t  eLT(uint, uint, uint64_t=0) override {
        const StateSimulator::label_1q_t L_I = StateSimulator::label_1q_t::I;
        return std::make_pair(L_I, L_I);
    }

    void snapshot(void) override {
    }

    void rollback_where(stim::simd_bits_range_ref) override {
    }
    
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
    cytnx::UniTensor qubit_base(int x) {
        cytnx::Tensor base({1, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = 1;

        cytnx::UniTensor out(base);
        out.set_labels({std::to_string(-(x+1)), std::to_string(x+1)});
        return out;
    }

    cytnx::Tensor id() {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        for (uint i = 0; i < W; i++) base.at({i, i}) = 1;
        return base;
    }

    cytnx::UniTensor h_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = sqrt(0.5);
        base.at({0, 1}) = sqrt(0.5);
        base.at({1, 0}) = sqrt(0.5);
        base.at({1, 1}) = -sqrt(0.5);
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor x_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 1}) = 1;
        base.at({1, 0}) = 1;
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor y_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 1}) = -M_J;
        base.at({1, 0}) = M_J;
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor z_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = 1;
        base.at({1, 1}) = -1;
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor s_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = 1;
        base.at({1, 1}) = M_J;
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor t_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = 1;
        base.at({1, 1}) = c64(cos(M_PI*0.25), sin(M_PI*0.25));
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor cx_base(uint x, uint y) {
        fp_t fault_angle;
        if (params.fix_faulty_cx_angle) {
            fault_angle = params.faulty_cx_angle;
        } else {
            fault_angle = rpdist(rng) * 2 * M_PI;
        }
        cytnx::Tensor base({W, W, W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0, 0, 0}) = 1;
        base.at({0, 1, 0, 1}) = 1;
        base.at({1, 0, 1, 1}) = 1;
        base.at({1, 1, 1, 0}) = 1;
        for (uint i = 2; i < W; i++) {
            // Apply RX error if either qubit is leaked.
            base.at({i, 0, i, 0}) = cos(fault_angle * 0.5);
            base.at({i, 0, i, 1}) = c64(0, -sin(fault_angle*0.5));
            base.at({i, 1, i, 0}) = c64(0, -sin(fault_angle*0.5));
            base.at({i, 1, i, 1}) = cos(fault_angle*0.5);

            base.at({0, i, 0, i}) = cos(fault_angle * 0.5);
            base.at({0, i, 1, i}) = c64(0, -sin(fault_angle*0.5));
            base.at({1, i, 0, i}) = c64(0, -sin(fault_angle*0.5));
            base.at({1, i, 1, i}) = cos(fault_angle*0.5);
        }
        return build_2q_gate(base, x, y);
    }

    cytnx::UniTensor rx_base(fp_t angle, uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = cos(angle*0.5);
        base.at({0, 1}) = c64(0, -sin(angle*0.5));
        base.at({1, 0}) = c64(0, -sin(angle*0.5));
        base.at({1, 1}) = cos(angle*0.5);
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor ry_base(fp_t angle, uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = cos(angle*0.5);
        base.at({0, 1}) = -sin(angle*0.5);
        base.at({1, 0}) = sin(angle*0.5);
        base.at({1, 1}) = cos(angle*0.5);
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor rz_base(fp_t angle, uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({0, 0}) = c64(cos(angle*0.5), -sin(angle*0.5));
        base.at({1, 1}) = c64(cos(angle*0.5), sin(angle*0.5));
        for (uint i = 2; i < W; i++) {
            base.at({i, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor m_base(uint x, uint b) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        base.at({b, b}) = 1;
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor reset_base(uint x) {
        cytnx::Tensor base({W, W}, cytnx::Type.ComplexDouble);
        for (uint i = 0; i < W; i++) {
            base.at({0, i}) = 1;
        }
        return build_1q_gate(base, x);
    }

    cytnx::UniTensor build_1q_gate(const cytnx::Tensor& base, uint x) {
        cytnx::UniTensor gate(base);
        // We must attach the entry and exit labels
        // for the tensor.
        std::string entry = qubit_to_label[x]; // This is the entry.
        std::string exit = std::to_string(label_ctr++);
        gate.set_labels({entry, exit});
        // Mark the exit label as the awaiting label for the qubit.
        awaiting_labels[x] = exit;
        return gate;
    }

    cytnx::UniTensor build_2q_gate(const cytnx::Tensor& base, uint x, uint y) {
        cytnx::UniTensor gate(base);
        // Attach entry and exit labels.
        std::string xentry = qubit_to_label[x];
        std::string yentry = qubit_to_label[y];
        std::string xexit = std::to_string(label_ctr++);
        std::string yexit = std::to_string(label_ctr++);
        gate.set_labels({xentry, yentry, xexit, yexit});
        // Mark the exit labels as awaiting for x and y.
        awaiting_labels[x] = xexit;
        awaiting_labels[y] = yexit;
        return gate;
    }

    cytnx::UniTensor test_1q_gate(const cytnx::UniTensor& gate, uint x) {
        std::string xlabel = qubit_to_label[x];
        cytnx::UniTensor* curr = label_to_tensor[xlabel];
        return Contract(*curr, gate);
    }

    cytnx::UniTensor test_2q_gate(const cytnx::UniTensor& gate, uint x, uint y) {
        std::string xlabel = qubit_to_label[x];
        std::string ylabel = qubit_to_label[y];
        cytnx::UniTensor* t1 = label_to_tensor[xlabel];
        cytnx::UniTensor* t2 = label_to_tensor[ylabel];

        cytnx::UniTensor s1 = Contract(*t1, gate);
        cytnx::UniTensor s2 = Contract(*t2, s1);
        return s2;
    }

    void apply_1q_gate(const cytnx::UniTensor& gate, uint x) {
        cytnx::UniTensor state = test_1q_gate(gate, x);
        update_1q_state(state, x);
    }

    void update_1q_state(const cytnx::UniTensor& state, uint x) {
        std::string xlabel = qubit_to_label[x];
        cytnx::UniTensor* curr = label_to_tensor[xlabel];
        cytnx::UniTensor* out = new cytnx::UniTensor(state);
        // Update data structures and delete curr.
        qubit_to_label[x] = awaiting_labels[x];
        label_to_qubit[awaiting_labels[x]] = x;
        label_to_qubit.erase(xlabel);
        label_to_tensor.erase(xlabel);

        tensor_pool.erase(curr);
        tensor_pool.insert(out);
        for (auto label : tensor_to_labels[curr]) {
            std::string true_label = label == xlabel ? awaiting_labels[x] : label;
            label_to_tensor[true_label] = out;
            tensor_to_labels[out].push_back(true_label);
        }
        tensor_to_labels.erase(curr);
        delete curr;
    }

    void apply_2q_gate(const cytnx::UniTensor& gate, uint x, uint y) {
        cytnx::UniTensor state = test_2q_gate(gate, x, y);
        update_2q_state(state, x, y);
    }

    void update_2q_state(const cytnx::UniTensor& state, uint x, uint y) {
        std::string xlabel = qubit_to_label[x];
        std::string ylabel = qubit_to_label[y];

        cytnx::UniTensor* t1 = label_to_tensor[xlabel];
        cytnx::UniTensor* t2 = label_to_tensor[ylabel];

        cytnx::UniTensor* out = new cytnx::UniTensor(state);

        // Update data structures and delete t1, t2.
        qubit_to_label[x] = awaiting_labels[x];
        qubit_to_label[y] = awaiting_labels[y];
        label_to_qubit[awaiting_labels[x]] = x;
        label_to_qubit[awaiting_labels[y]] = y;
        label_to_qubit.erase(xlabel);
        label_to_qubit.erase(ylabel);

        label_to_tensor.erase(xlabel);
        label_to_tensor.erase(ylabel);

        tensor_pool.erase(t1);
        tensor_pool.erase(t2);
        tensor_pool.insert(out);
        std::vector<std::string> t1_t2_labels;
        if (t1 == t2) {
            t1_t2_labels = tensor_to_labels[t1];
            tensor_to_labels.erase(t1);
            delete t1;
        } else {
            t1_t2_labels = tensor_to_labels[t1];
            for (auto s : tensor_to_labels[t2]) {
                t1_t2_labels.push_back(s);
            }
            tensor_to_labels.erase(t1);
            tensor_to_labels.erase(t2);
            delete t1;
            delete t2;
        }
        for (auto label : t1_t2_labels) {
            std::string true_label;
            if (label == xlabel)        true_label = awaiting_labels[x];
            else if (label == ylabel)   true_label = awaiting_labels[y];
            else {
                // Apply identity.
                cytnx::UniTensor u(id());
                true_label = std::to_string(label_ctr++);
                u.set_labels({label, true_label});
                *out = Contract(*out, u);
                // Update data structures.
                label_to_tensor.erase(label);

                label_to_qubit[true_label] = label_to_qubit[label];
                qubit_to_label[label_to_qubit[label]] = true_label;
                label_to_qubit.erase(label);
            }
            label_to_tensor[true_label] = out;
            tensor_to_labels[out].push_back(true_label);
        }
        // Check if the bond size is too high.
        if (is_bond_length_critical())  compress();
    }

    bool is_bond_length_critical(void) { return false; }
    void compress(void) {}

    std::map<uint, std::string>                 qubit_to_label;
    std::map<std::string, uint>                 label_to_qubit;
    std::set<cytnx::UniTensor*>                 tensor_pool;
    std::map<std::string, cytnx::UniTensor*>    label_to_tensor;
    std::map<cytnx::UniTensor*, std::vector<std::string>>   tensor_to_labels;

    std::map<uint, std::string> awaiting_labels;    // New labels after operation.

    int32_t label_ctr;

    std::uniform_real_distribution<fp_t>    rpdist;
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

template <int b> cytnx::Tensor
build_rotation_operator_for_b2_subspace(fp_t, fp_t);

}   // qontra

#endif  // UNIVERSAL_SIM_h
