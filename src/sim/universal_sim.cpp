/*
 *  author: Suhas Vittal
 *  date:   21 August 2023
 * */

#include "sim/universal_sim.h"

using namespace cytnx;

typedef std::complex<double> c64;

namespace qontra {

template <int W>
UniversalSimulator::UniversalSimulator(uint n)
    :StateSimulator(n, 1),
    qubit_to_label(),
    tensor_pool(),
    label_to_tensor(),
    tensor_to_labels(),
    awaiting_labels(),
    label_ctr(n+1)
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

template <int W> UniTensor
UniversalSimulator::h_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({0, 1}) = 1;
    base.at({1, 0}) = 1;
    base.at({1, 1}) = -1;
    base = base * sqrt(0.5);
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::x_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 1}) = 1;
    base.at({1, 0}) = 1;
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::y_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 1}) = c64(0, -1);
    base.at({1, 0}) = c64(0, 1);
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::z_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = -1;
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::s_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = c64(0, 1);
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::t_base(uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = 1;
    base.at({1, 1}) = c64(cos(M_PI*0.25), sin(M_PI*0.25));
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::cx_base(uint x, uint y) {
    Tensor base({2, 2, 2, 2}, type=Type.ComplexDouble);
    base.at({0, 0, 0, 0}) = 1;
    base.at({0, 1, 0, 1}) = 1;
    base.at({1, 0, 1, 1}) = 1;
    base.at({1, 1, 1, 0}) = 1;
    return build_2q_gate(base, x, y);
}

template <int W> UniTensor
UniversalSimulator::rx_base(fp_t angle, uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = cos(angle*0.5);
    base.at({0, 1}) = c64(0, -sin(angle*0.5));
    base.at({1, 0}) = c64(0, -sin(angle*0.5));
    base.at({1, 1}) = cos(angle*0.5);
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::ry_base(fp_t angle, uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = cos(angle*0.5);
    base.at({0, 1}) = -sin(angle*0.5);
    base.at({1, 0}) = sin(angle*0.5);
    base.at({1, 1}) = cos(angle*0.5);
    return build_1q_gate(base, x);
}

template <int W> UniTensor
UniversalSimulator::rz_base(fp_t angle, uint x) {
    Tensor base({2, 2}, type=Type.ComplexDouble);
    base.at({0, 0}) = c64(cos(angle*0.5), -sin(angle*0.5));
    base.at({1, 1}) =  c64(cos(angle*0.5), sin(angle*0.5));
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

template <int W> void
UniversalSimulator::apply_1q_gate(const UniTensor gate, uint x) {
    std::string xlabel = qubit_to_label[x];
    UniTensor* curr = label_to_tensor[xlabel];

    UniTensor* outcome = new UniTensor(Contract(curr, gate));

    // Update data structures and delete curr.
    qubit_to_label[x] = awaiting_labels[x];
    label_to_tensor.erase(xlabel);

    tensor_pool.erase(curr);
    tensor_pool.insert(outcome);
    for (auto label : tensor_to_labels[curr]) {
        auto true_label = label == xlabel ? awaiting_labels[x] : label;
        label_to_tensor[true_label] = outcome;
        tensor_to_label[outcome].push_back(true_label);
    }
    tensor_to_label.erase(curr);
    delete curr;
    // Note that this operation cannot increase bond length.
}

template <int W> void
UniversalSimulator::apply_2q_gate(const UniTensor gate, uint x, uint y) {
    std::string xlabel = qubit_to_label[x];
    std::string ylabel = qubit_to_label[y];
    UniTensor* t1 = label_to_tensor(xlabel);
    UniTensor* t2 = label_to_tensor(ylabel);

    UniTensor s1 = Contract(t1, gate);
    UniTensor s2 = Contract(t2, s1);
    UniTensor* outcome = new UniTensor(s2);

    // Update data structures and delete t1, t2.
    qubit_to_label[x] = awaiting_labels[x];
    qubit_to_label[y] = awaiting_labels[y];
    label_to_tensor.erase(x);
    label_to_tensor.erase(y);

    tensor_pool.erase(t1);
    tensor.pool.erase(t2);
    tensor_pool.insert(outcome);
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
        label_to_tensor[true_label] = outcome;
        tensor_to_labels[outcome].push_back(true_label);
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

}   // qontra
