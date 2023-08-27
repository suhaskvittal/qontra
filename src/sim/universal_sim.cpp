/*
 *  author: Suhas Vittal
 *  date:   21 August 2023
 * */

#include "sim/universal_sim.h"

using namespace cytnx;

typedef std::complex<double> c64;

namespace qontra {

#define L_I     StateSimulator::label_1q_t::I
#define L_X     StateSimulator::label_1q_t::X
#define L_Y     StateSimulator::label_1q_t::Y
#define L_Z     StateSimulator::label_1q_t::Z
#define L_L     StateSimulator::label_1q_t::L

static const StateSimulator::label_1q_t IXYZ[4] = { L_I, L_X, L_Y, L_Z };
static const StateSimulator::label_1q_t IX[2] = { L_I, L_X };
static const StateSimulator::label_1q_t IY[2] = { L_I, L_Y };
static const StateSimulator::label_1q_t IZ[2] = { L_I, L_Z };
static const StateSimulator::label_1q_t IL[2] = { L_I, L_L };

static const c64 M_J = c64(0, 1);

template <int W>
UniversalSimulator::UniversalSimulator(uint n)
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
    for (int i = 0; i < n; i++) {
        std::string label = std::to_string(i+1);
        UniTensor* init = new UniTensor(qubit_base(i));

        qubit_to_label[i] = label;
        label_to_tensor[label] = init;
        tensor_to_labels[init].push_back(label);
        tensor_pool.insert(init);
    }
}

template <int W> void
UniversalSimulator::M(std::vector<uint> operands, fp_t m1w0, fp_t m0w1, int record) {
    for (uint x : operands) {
        fp_t r = rpdist(rng);
        // Get the probability of measuring each state.
        for (uint i = 0; i < W; i++) {
            auto tmp = test_1q_gate(m_base(x, i), x);
            fp_t pr = (fp_t) Scalar(tmp.Norm().at({0}).abs());
            if (r < pr) {
                // Collapse to this state.
                tmp.normalize_();
                update_1q_state(tmp);
                if (record >= 0) {
                    uint m = i;
                    // Randomly resolve the measurement if we collapsed to |2>, |3>, etc.
                    if (m > 1) {
                        fp_t x = rpdist(rng);
                        m = x < 0.5 ? 0 : 1;
                    }
                    // Apply errors.
                    fp_t e = rpdist(rng);
                    if (m == 0 && e < m1w0)         m = 1;
                    else if (m == 1 && e < m0w1)    m = 0;

                    record_table[record++] = m;
                }
                break;
            } else {
                r -= pr;
            }
        }
    }
}

template <int W> StateSimulator::label_1q_t
UniversalSimulator::eDP1(uint x, uint64_t t) {
    auto p = rng() & 3;
    if (p & 1)  apply_1q_gate(x_base(x), x);
    if (p & 2)  apply_2q_gate(z_base(x), x);
    return IXYZ[p];
}

template <int W> StateSimulator::label_1q_t
UniversalSimulator::eX(uint x, uint64_t t) {
    apply_1q_gate(x_base(x), x);
    return L_X;
}

template <int W> StateSimulator::label_1q_t
UniversalSimulator::eY(uint x, uint64_t t) {
    apply_1q_gate(y_base(x), x);
    return L_Y;
}

template <int W> StateSimulator::label_1q_t
UniversalSimulator::eZ(uint x, uint64_t t) {
    apply_1q_gate(z_base(x), x);
    return L_Z;
}

template <int W> StateSimulator::label_1q_t
UniversalSimulator::eL(uint x, uint64_t t) {
    // Leave this for explicitly specialized implementations.
    return L_I;
}

template <int W, int b1, int b2> StateSimulator::label_1q_t
UniversalSimulator::eAD(uint x, uint64_t t) {
    fp_t angle;
    if (params.fix_amplitude_damping_angle) {
        angle = params.amplitude_damping_angle;
    } else {
        angle = rpdist(rng) * 2 * M_PI;
    }

    Tensor rot({W, W}, type=Type.ComplexDouble);
    for (int i = 0; i < W; i++) {
        if (i != b1 && i != b2) rot.at({i, i}) = 1;
    }
    rot.at({b1, b1}) = 1.0;
    rot.at({b1, b2}) = c64(0.0, -sin(angle*0.5));
    rot.at({b2, b2}) = cos(angle*0.5);

    apply_1q_gate(build_1q_gate(rot, x), x);
    return b1 <= 1 && b2 <= 1 ? L_X : L_L;
}

template <int W, int b1, int b2> StateSimulator::label_1q_t
UniversalSimulator::ePD(uint x, uint64_t t) {
    fp_t angle;
    if (params.fix_phase_damping_angle) {
        angle = params.phase_damping_angle;
    } else {
        angle = rpdist(rng) * 2 * M_PI;
    }

    Tensor rot({W, W}, type=Type.ComplexDouble);
    for (int i = 0; i < W; i++) {
        if (i != b1 && i != b2) rot.at({i, i}) = 1;
    }
    rot.at({b1, b1}) = c64(cos(angle*0.5), -sin(angle*0.5));
    rot.at({b2, b2}) = c64(cos(angle*0.5), sin(angle*0.5));
    
    apply_1q_gate(build_1q_gate(rot, x), x);
    return L_Z;
}

template <int W, int b1, int b2> StateSimulator::label_1q_t
UniversalSimulator::eHEAT(uint x, uint64_t t) {
    fp_t angle;
    if (params.fix_heating_angle) {
        angle = params.heating_angle;
    } else {
        angle = rpdist(rng) * 2 * M_PI;
    }

    Tensor rot({W, W}, type=Type.ComplexDouble);
    for (int i = 0; i < W; i++) {
        if (i != b1 && i != b2) rot.at({i, i}) = 1;
    }
    rot.at({b1, b1}) = cos(angle*0.5);
    rot.at({b2, b1}) = c64(0.0, -sin(angle*0.5));
    rot.at({b2, b2}) = 1.0;

    apply_1q_gate(build_1q_gate(rot, x), x);
    return b1 <= 1 && b2 <= 1 ? L_X : L_L;
}

template <int W> StateSimulator::label_2q_t
UniversalSimulator::eDP2(uint x, uint y, uint64_t t) {
    auto p = rng() & 15;
    if (p & 1)  apply_1q_gate(x_base(x), x);
    if (p & 2)  apply_1q_gate(z_base(x), x);
    if (p & 4)  apply_1q_gate(x_base(y), y);
    if (p & 8)  apply_1q_gate(z_base(y), y);
    return std::make_pair(IXYZ[p & 3], IXYZ[(p>>2) & 3]);
}

template <int W> StateSimulator::label_2q_t
UniversalSimulator::eLI(uint x, uint y, uint64_t t) {
    // Again, as with eL, leave this to be specialized.
    return L_I;
}

template <int W> StateSimulator::label_2q_t
UniversalSimulator::eLT(uint x, uint y, uint64_t t) {
    // Again, as with eL, leave this to be specialized.
    return L_I;
}

template <int W> UniTensor
UniversalSimulator::qubit_base(int x) {
    Tensor base({1, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;

    UniTensor out(base);
    out.set_labels({std::to_string(-(x+1)), std::to_string(x+1)});
    return out;
}

template <int W> UniTensor
UniversalSimulator::h_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = sqrt(0.5);
    base.at({0, 1}) = sqrt(0.5);
    base.at({1, 0}) = sqrt(0.5);
    base.at({1, 1}) = -sqrt(0.5);
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::x_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 1}) = 1;
    base.at({1, 0}) = 1;
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::y_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 1}) = -M_J;
    base.at({1, 0}) = M_J;
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::z_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = -1;
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::s_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = M_J;
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::t_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = c64(cos(M_PI*0.25), sin(M_PI*0.25));
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::cx_base(uint x, uint y) {
    fp_t fault_angle;
    if (params.fix_faulty_cx_angle) {
        fault_angle = params.faulty_cx_angle;
    } else {
        fault_angle = rpdist(rng) * 2 * M_PI;
    }
    Tensor base({W, W, W, W}, type=Type.ComplexDouble);
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

template <int W> UniTensor
UniversalSimulator::rx_base(fp_t angle, uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = cos(angle*0.5);
    base.at({0, 1}) = c64(0, -sin(angle*0.5));
    base.at({1, 0}) = c64(0, -sin(angle*0.5));
    base.at({1, 1}) = cos(angle*0.5);
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::ry_base(fp_t angle, uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = cos(angle*0.5);
    base.at({0, 1}) = -sin(angle*0.5);
    base.at({1, 0}) = sin(angle*0.5);
    base.at({1, 1}) = cos(angle*0.5);
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::rz_base(fp_t angle, uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({0, 0}) = c64(cos(angle*0.5), -sin(angle*0.5));
    base.at({1, 1}) =  c64(cos(angle*0.5), sin(angle*0.5));
    for (uint i = 2; i < W; i++) {
        base.at({i, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::m_base(uint x, uint b) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    base.at({b, b}) = 1;
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::reset_base(uint x) {
    Tensor base({W, W}, type=Type.ComplexDouble);
    for (uint i = 0; i < W; i++) {
        base.at({0, i}) = 1;
    }
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::build_1q_gate(const Tensor& base, uint x) {
    UniTensor gate(base);
    // We must attach the entry and exit labels
    // for the tensor.
    std::string entry = qubit_to_label[x]; // This is the entry.
    std::string exit = std::to_string(label_ctr++);
    gate.set_labels({entry, exit});
    // Mark the exit label as the awaiting label for the qubit.
    awaiting_labels[x] = exit;
    return gate;
}

template <int W> UniTensor
UniversalSimulator::build_2q_gate(const Tensor& base, uint x, uint y) {
    UniTensor gate(base);
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

template <int W> UniTensor
UniversalSimulator::test_1q_gate(const cytnx::UniTensor gate, uint x) {
    std::string xlabel = qubit_to_label[x];
    UniTensor* curr = label_to_tensor[xlabel];
    return Contract(curr, gate);
}

template <int W> UniTensor
UniversalSimulator::test_2q_gate(const cytnx::UniTensor gate, uint x, uint y) {
    std::string xlabel = qubit_to_label[x];
    std::string ylabel = qubit_to_label[y];
    UniTensor* t1 = label_to_tensor(xlabel);
    UniTensor* t2 = label_to_tensor(ylabel);

    UniTensor s1 = Contract(t1, gate);
    UniTensor s2 = Contract(t2, s1);
    return s2;
}

template <int W> void
UniversalSimulator::apply_1q_gate(const UniTensor gate, uint x) {
    update_1q_state(test_1q_gate(gate, x), x);
}

template <int W> void
UniversalSimulator::update_1q_state(const UniTensor state, uint x) {
    std::string xlabel = qubit_to_label[x];

    UniTensor* out = new UniTensor(state);
    // Update data structures and delete curr.
    qubit_to_label[x] = awaiting_labels[x];
    label_to_tensor.erase(xlabel);

    tensor_pool.erase(curr);
    tensor_pool.insert(out);
    for (auto label : tensor_to_labels[curr]) {
        auto true_label = label == xlabel ? awaiting_labels[x] : label;
        label_to_tensor[true_label] = outcome;
        tensor_to_label[out].push_back(true_label);
    }
    tensor_to_label.erase(curr);
    delete curr;
}

template <int W> void
UniversalSimulator::apply_2q_gate(const UniTensor gate, uint x, uint y) {
    update_2q_state(test_2q_gate(gate, x, y), x, y);
}

template <int W> void
UniversalSimulator::update_2q_state(const UniTensor state, uint x, uint y) {
    std::string xlabel = qubit_to_label[x];
    std::string ylabel = qubit_to_label[y];

    UniTensor* out = new UniTensor(state);

    // Update data structures and delete t1, t2.
    qubit_to_label[x] = awaiting_labels[x];
    qubit_to_label[y] = awaiting_labels[y];
    label_to_tensor.erase(x);
    label_to_tensor.erase(y);

    tensor_pool.erase(t1);
    tensor.pool.erase(t2);
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
        else                        true_label = label;
        label_to_tensor[true_label] = out;
        tensor_to_labels[out].push_back(true_label);
    }
    // Check if the bond size is too high.
    if (is_bond_length_critical())  compress();
}

// TODO:

template <int W> bool
UniversalSimulator::is_bond_length_critical() {
    return false;
}

template <int W> void
UniversalSimulator::compress() {
}

//
//  TEMPLATE SPECIALIZATIONS FOR W = 3
//

template <> StateSimulator::label_2q_t
UniversalSimulator<3>::eLI(uint x, uint y, uint64_t t) {
    apply_1q_gate(build_1q_gate(build_leakage_operator(), x), x);
    apply_1q_gate(build_1q_gate(build_leakage_operator(), y), y);
    return std::make_pair(L_L, L_L);
}

template <> StateSimulator::label_2q_t
UniversalSimulator<3>::eLT(uint x, uint y, uint64_t t) {
    fp_t transport_angle;
    if (params.fix_leakage_transport_angle) {
        transport_angle = params.leakage_transport_angle;
    } else {
        transport_angle = rpdist(rng) * 2 * M_PI;
    }
    // We will apply the transport across four subspaces.
    //  a. |02> - |22>
    //  b. |12> - |22>
    //  c. |20> - |22>
    //  d. |21> - |22>
    //
    typedef std::tuple<int, int>  subspace_t;
    const std::vector<subspace_t> transport_subspaces
    {
        std::make_tuple(0, 2),
        std::make_tuple(1, 2),
        std::make_tuple(2, 0),
        std::make_tuple(2, 1)
    };

    Tensor transport_operator({3, 3, 3, 3}, type=Type.ComplexDouble);
    for (auto ss : transport_subspaces) {
        int x = std::get<0>(ss), y = std::get<1>(ss);
        transport_operator.at({x, 2, x, 2}) = cos(transport_angle*0.5);
        transport_operator.at({x, 2, y, 2}) = -sin(transport_angle*0.5);
        transport_operator.at({y, 2, x, 2}) = sin(transport_angle*0.5);
        transport_operator.at({y, 2, y, 2}) = cos(transport_angle*0.5);
    }
    apply_2q_gate(build_2q_gate(transport_operator, x, y), x, y);
    // We can't accurately say what happened, so just return identity
    // labels.
    return std::make_pair(L_I, L_I);
}

//
// Helper functions
//

Tensor
build_leakage_operator() {
    fp_t a1 = rpdist(rng) * 2 * M_PI;
    fp_t a2 = rpdist(rng) * 2 * M_PI;
    fp_t a3 = rpdist(rng) * 2 * M_PI;
    
    Tensor r02 = build_rotation_operator_for_b2_subspace<0>(a1, a2);
    Tensor r12 = build_rotation_operator_for_b2_subspace<1>(a1, a2);
    Tensor uz({3, 3});
    uz.at({0, 0}) = 1;
    uz.at({1, 1}) = 1;
    uz.at({2, 2}) = c64(cos(a3), sin(a3));

    return uz * r02 * r12;
}

template <int b> Tensor
build_rotation_operator_for_b2_subspace(fp_t a1, fp_t a2) {
    using namespace linalg;

    Tensor pauli_x({3, 3}, type=Type.ComplexDouble);
    pauli_x.at({b, 2}) = 1;
    pauli_x.at({2, b}) = 1;
    pauli_x.at({1-b, 1-b}) = 1;

    Tensor pauli_y({3, 3}, type=Type.ComplexDouble);
    pauli_x.at({b, 2}) = -M_J;
    pauli_x.at({2, b}) = M_J;
    pauli_x.at({1-b, 1-b}) = 1;

    Tensor cos_a2X = 0.5*a2*(ExpH(M_J*pauli_x) + ExpH(-M_J*pauli_x));
    Tensor sin_a2Y = 0.5*M_J*a2*(ExpH(M_J*pauli_y) - Exp_H(-M_J*pauli_y));
    
    Tensor u = c64(cos(0.5*a1), sin(0.5*a2)) * ExpH(0.5*M_J*a1 * (cos_a2X+sin_a2Y));
    return u;
}

}   // qontra
