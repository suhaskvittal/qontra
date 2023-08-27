/*
 *  author: Suhas Vittal
 *  date:   21 August 2023
 * */

#include "sim/universal_sim.h"

using namespace cytnx;

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

    Tensor transport_operator({3, 3, 3, 3}, Type.ComplexDouble);
    for (auto ss : transport_subspaces) {
        uint64_t x = std::get<0>(ss), y = std::get<1>(ss);
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

static std::mt19937_64 __rng(0);
static std::uniform_real_distribution __rangdist(0.0, 2*M_PI);

Tensor
build_leakage_operator() {
    fp_t a1 = __rangdist(__rng);
    fp_t a2 = __rangdist(__rng);
    fp_t a3 = __rangdist(__rng);
    
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

    Tensor pauli_x({3, 3}, Type.ComplexDouble);
    pauli_x.at({b, 2}) = 1;
    pauli_x.at({2, b}) = 1;
    pauli_x.at({1-b, 1-b}) = 1;

    Tensor pauli_y({3, 3}, Type.ComplexDouble);
    pauli_x.at({b, 2}) = -M_J;
    pauli_x.at({2, b}) = M_J;
    pauli_x.at({1-b, 1-b}) = 1;

    Tensor cos_a2X = 0.5*a2*(ExpH(M_J*pauli_x) + ExpH(-M_J*pauli_x));
    Tensor sin_a2Y = 0.5*M_J*a2*(ExpH(M_J*pauli_y) - ExpH(-M_J*pauli_y));
    
    Tensor u = c64(cos(0.5*a1), sin(0.5*a2)) * ExpH(0.5*M_J*a1 * (cos_a2X+sin_a2Y));
    return u;
}

}   // qontra
