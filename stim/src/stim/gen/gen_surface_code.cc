#include "stim/gen/gen_surface_code.h"

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>

using namespace stim;

#define USE_SWAPS

struct surface_coord {
    float x;
    float y;
    surface_coord operator+(surface_coord other) const {
        return {x + other.x, y + other.y};
    }
    surface_coord operator-(surface_coord other) const {
        return {x - other.x, y - other.y};
    }
    bool operator==(surface_coord other) const {
        return x == other.x && y == other.y;
    }
    bool operator<(surface_coord other) const {
        if (x != other.x) {
            return x < other.x;
        }
        return y < other.y;
    }
};

GeneratedCircuit _finish_surface_code_circuit(
    std::function<uint32_t(surface_coord)> coord_to_index,
    const std::set<surface_coord> &data_coords,
    const std::set<surface_coord> &x_measure_coords,
    const std::set<surface_coord> &z_measure_coords,
    const CircuitGenParameters &params,
    const std::vector<surface_coord> &x_order,
    const std::vector<surface_coord> &z_order,
    const std::vector<surface_coord> x_observable,
    const std::vector<surface_coord> z_observable,
    bool is_memory_x) 
{
    if (params.rounds < 1) {
        throw std::invalid_argument("Need rounds >= 1.");
    }
    if (params.distance < 2) {
        throw std::invalid_argument("Need a distance >= 2.");
    }

    const auto &chosen_basis_observable = is_memory_x ? x_observable : z_observable;
    const auto &chosen_basis_measure_coords = is_memory_x ? x_measure_coords : z_measure_coords;

    // Index the measurement qubits and data qubits.
    std::map<surface_coord, uint32_t> p2q;
    for (auto q : data_coords) {
        p2q[q] = coord_to_index(q);
    }
    for (auto q : x_measure_coords) {
        p2q[q] = coord_to_index(q);
    }
    for (auto q : z_measure_coords) {
        p2q[q] = coord_to_index(q);
    }

    // Reverse index.
    std::map<uint32_t, surface_coord> q2p;
    for (const auto &kv : p2q) {
        q2p[kv.second] = kv.first;
    }

    // Make target lists for various types of qubits.
    std::vector<uint32_t> data_qubits;
    std::vector<uint32_t> measurement_qubits;
    std::vector<uint32_t> x_measurement_qubits;
    std::vector<uint32_t> z_measurement_qubits;
    std::vector<uint32_t> all_qubits;
    for (auto q : data_coords) {
        data_qubits.push_back(p2q[q]);
    }
    for (auto q : x_measure_coords) {
//        measurement_qubits.push_back(p2q[q]);
        x_measurement_qubits.push_back(p2q[q]);
    }
    for (auto q : z_measure_coords) {
//        measurement_qubits.push_back(p2q[q]);
        z_measurement_qubits.push_back(p2q[q]);
    }
    all_qubits.insert(all_qubits.end(), data_qubits.begin(), data_qubits.end());
    all_qubits.insert(all_qubits.end(), measurement_qubits.begin(), measurement_qubits.end());
    std::sort(all_qubits.begin(), all_qubits.end());
    std::sort(data_qubits.begin(), data_qubits.end());
//    std::sort(measurement_qubits.begin(), measurement_qubits.end());
    std::sort(x_measurement_qubits.begin(), x_measurement_qubits.end());
    std::sort(z_measurement_qubits.begin(), z_measurement_qubits.end());
    for (uint32_t xm : x_measurement_qubits) {
        measurement_qubits.push_back(xm);
    }
    for (uint32_t zm : z_measurement_qubits) {
        measurement_qubits.push_back(zm);
    }

    // Reverse index the measurement order used for defining detectors.
    std::map<surface_coord, uint32_t> data_coord_to_order;
    std::map<surface_coord, uint32_t> measure_coord_to_order;
    for (auto q : data_qubits) {
        auto i = data_coord_to_order.size();
        data_coord_to_order[q2p[q]] = i;
    }
    for (auto q : measurement_qubits) {
        auto i = measure_coord_to_order.size();
        measure_coord_to_order[q2p[q]] = i;
    }

    // List out CNOT gate targets using given interaction orders.
    std::array<std::vector<uint32_t>, 4> cnot_targets;
    // Only use this if using SWAP LRU.
    std::array<std::vector<uint32_t>, 2> first_round_swap_targets;
    // If we are using SWAP LRU, then every 4 rounds of syndrome
    // extraction are different. Furthermore, the first round is
    // different from the remaining rounds.
    // 
    // We only perform SWAPs on Z stabilizers because the overhead
    // is 1 CNOT. In contrast, X stabilizers require a whole SWAP
    // (3 CNOTs).
    const uint32_t INVALID_SWAP = (uint32_t)-1;
    const uint32_t lru_cycles = params.swap_lru_with_no_swap ? 5 : 4;
    std::map<uint32_t, std::array<uint32_t, 5>> swap_order;
    if (params.use_swap_lru) {
        // Build swap order.
        for (size_t cycle = 0; cycle < 4; cycle++) {
            std::set<uint32_t> already_swapped;
            for (auto measure : z_measure_coords) {
                uint32_t pm = p2q[measure];
                if (!swap_order.count(pm)) {
                    swap_order[pm] = std::array<uint32_t, 5>();
                }

                bool found = true;
                size_t i = cycle;
                auto data = measure + z_order[i];
                i = (i+1) & 0x3;
                while (p2q.find(data) == p2q.end() 
                        || already_swapped.count(p2q[data]))
                {
                    data = measure + z_order[i];
                    i = (i + 1) & 0x3;
                    if (i == cycle) {
                        // Then we have been unable to find a proper
                        // swap candidate. This ancilla will not be
                        // swapping for this cycle.
                        found = false;
                        break;
                    }
                }

                if (found) {
                    swap_order[pm][cycle] = p2q[data];
                    already_swapped.insert(p2q[data]);
                } else {
                    swap_order[pm][cycle] = INVALID_SWAP;
                }
            }
        }

        if (params.swap_lru_with_no_swap) {
            for (auto measure : z_measure_coords) {
                auto pm = p2q[measure];
                swap_order[pm][4] = pm;
            }
        }
    }

    for (size_t k = 0; k < 4; k++) {
        for (auto measure : x_measure_coords) {
            auto data = measure + x_order[k];
            if (p2q.find(data) != p2q.end()) {
                cnot_targets[k].push_back(p2q[measure]);
                cnot_targets[k].push_back(p2q[data]);
            }
        }

        for (auto measure : z_measure_coords) {
            uint32_t pm = p2q[measure];
            auto data = measure + z_order[k];
            if (p2q.find(data) != p2q.end()) {
                uint32_t pd = p2q[data];
                if (params.use_swap_lru && pd == swap_order[pm][0]) {
                    first_round_swap_targets[0].push_back(pm);
                    first_round_swap_targets[0].push_back(swap_order[pm][0]);
                    first_round_swap_targets[1].push_back(swap_order[pm][0]);
                    first_round_swap_targets[1].push_back(pm);
#ifndef USE_SWAPS
                    continue;
#endif
                }
                cnot_targets[k].push_back(pd);
                cnot_targets[k].push_back(pm);
            }
        }
    }

    // Build the repeated actions that make up the surface code cycle.
    Circuit cycle_actions;  // Use as the first round of syndrome
                            // extraction if using SWAP LRU.
    params.append_begin_round_tick(cycle_actions, data_qubits);
    params.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    for (const auto &targets : cnot_targets) {
        cycle_actions.append_op("TICK", {});
        params.append_unitary_2(cycle_actions, "CNOT", targets);
    }
    // Introduce SWAP operations into surface code cycle.
    if (params.use_swap_lru) {
        cycle_actions.append_op("TICK", {});
#ifdef USE_SWAPS
        params.append_unitary_2(cycle_actions, "SWAP", first_round_swap_targets[0]);
#else
        params.append_unitary_2(cycle_actions, "CNOT", first_round_swap_targets[0]);
        cycle_actions.append_op("TICK", {});
        params.append_unitary_2(cycle_actions, "CNOT", first_round_swap_targets[1]);
#endif
    }
    cycle_actions.append_op("TICK", {});
    params.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    cycle_actions.append_op("TICK", {});
    if (params.use_swap_lru) {
        params.append_measure_reset(cycle_actions, x_measurement_qubits);
        std::vector<uint32_t> local_measurement_qubits;
        for (uint32_t m : z_measurement_qubits) {
            if (swap_order[m][0] != INVALID_SWAP) {
                local_measurement_qubits.push_back(swap_order[m][0]);
            } else {
                local_measurement_qubits.push_back(m);
            }
        }
        cycle_actions.append_op("TICK", {});
        params.append_measure_reset(cycle_actions, local_measurement_qubits);
    } else {
        params.append_measure_reset(cycle_actions, measurement_qubits);
    }

    std::array<Circuit, 5> alt_cycle_actions;
    if (params.use_swap_lru) {
        for (size_t cycle = 0; cycle < lru_cycles; cycle++) {
            const size_t prev_cycle = cycle == 0 
                                         ? lru_cycles-1 : cycle - 1;

            std::array<std::vector<uint32_t>, 4> alt_cnot_targets;
            std::vector<uint32_t> alt_pre_swap_targets;
            std::array<std::vector<uint32_t>, 2> alt_post_swap_targets;
            for (size_t k = 0; k < 4; k++) {
                for (auto measure : x_measure_coords) {
                    auto data = measure + x_order[k];
                    if (p2q.find(data) != p2q.end()) {
                        alt_cnot_targets[k].push_back(p2q[measure]);
                        alt_cnot_targets[k].push_back(p2q[data]);
                    }
                }

                for (auto measure : z_measure_coords) {
                    uint32_t pm = p2q[measure];
                    auto data = measure + z_order[k];
                    if (p2q.find(data) != p2q.end()) {
                        uint32_t pd = p2q[data];
                        if (pd == swap_order[pm][prev_cycle]) {
                            alt_pre_swap_targets.push_back(pm);
                            alt_pre_swap_targets.push_back(swap_order[pm][prev_cycle]);
#ifndef USE_SWAPS
                            continue;
#endif
                        }
                        if (pd == swap_order[pm][cycle]) {
                            alt_post_swap_targets[0].push_back(pm);
                            alt_post_swap_targets[0].push_back(swap_order[pm][cycle]);
                            alt_post_swap_targets[1].push_back(swap_order[pm][cycle]);
                            alt_post_swap_targets[1].push_back(pm);
#ifndef USE_SWAPS
                            continue;
#endif 
                        }
                        alt_cnot_targets[k].push_back(p2q[data]);
                        alt_cnot_targets[k].push_back(pm);
                    }
                }
            }
            Circuit alt_cycle;
            params.append_begin_round_tick(alt_cycle, data_qubits);
            params.append_unitary_1(alt_cycle, "H", x_measurement_qubits);
            if (!alt_pre_swap_targets.empty()) {
#ifdef USE_SWAPS
                params.append_unitary_2(alt_cycle, "SWAP", alt_pre_swap_targets);
#else
                params.append_unitary_2(alt_cycle, "CNOT", alt_pre_swap_targets);
#endif
            }
            for (const auto &targets : cnot_targets) {
                alt_cycle.append_op("TICK", {});
                params.append_unitary_2(alt_cycle, "CNOT", targets);
            }
            alt_cycle.append_op("TICK", {});
            params.append_unitary_1(alt_cycle, "H", x_measurement_qubits);
            if (!alt_post_swap_targets[0].empty()) {
                alt_cycle.append_op("TICK", {});
#ifdef USE_SWAPS
                params.append_unitary_2(alt_cycle, "SWAP", alt_post_swap_targets[0]);
#else
                params.append_unitary_2(alt_cycle, "CNOT", alt_post_swap_targets[0]);
                alt_cycle.append_op("TICK", {});
                params.append_unitary_2(alt_cycle, "CNOT", alt_post_swap_targets[1]);
#endif
            }
            alt_cycle.append_op("TICK", {});
            params.append_measure_reset(alt_cycle, x_measurement_qubits);
            std::vector<uint32_t> local_measurement_qubits;
            for (uint32_t m : z_measurement_qubits) {
                if (swap_order[m][cycle] != INVALID_SWAP) {
                    local_measurement_qubits.push_back(swap_order[m][cycle]);
                } else {
                    local_measurement_qubits.push_back(m);
                }
            }
            alt_cycle.append_op("TICK", {});
            params.append_measure_reset(alt_cycle, local_measurement_qubits);

            alt_cycle_actions[cycle] = alt_cycle;
        }
    }

    // Build the start of the circuit, getting a state that's ready to cycle.
    // In particular, the first cycle has different detectors and so has to be handled special.
    Circuit head;
    for (const auto &kv : q2p) {
        head.append_op("QUBIT_COORDS", {kv.first}, {kv.second.x, kv.second.y});
    }
    params.append_reset(head, data_qubits, "ZX"[is_memory_x]);
    params.append_reset(head, measurement_qubits);
    head += cycle_actions;
    for (auto measure : chosen_basis_measure_coords) {
        head.append_op(
            "DETECTOR",
            {(uint32_t)(measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT},
            {measure.x, measure.y, 0});
    }

    // Build the repeated body of the circuit, including the detectors comparing to previous cycles.
    std::array<Circuit, 5> circuit_bodies;
    if (params.use_swap_lru) {
        for (size_t cycle = 0; cycle < lru_cycles; cycle++) {
            circuit_bodies[cycle] = alt_cycle_actions[cycle];
        }
    } else {
        circuit_bodies[0] = cycle_actions;
    }
    for (size_t cycle = 0; cycle < lru_cycles; cycle++) {
        if (!params.use_swap_lru && cycle > 0) {
            break;
        }
        Circuit& body = circuit_bodies[cycle];
        body.append_op("SHIFT_COORDS", {}, {0, 0, 1});
        uint32_t m = measurement_qubits.size();
        if (params.both_stabilizers) {
            for (auto m_index : measurement_qubits) {
                auto m_coord = q2p[m_index];
                auto k = (uint32_t)measurement_qubits.size() - measure_coord_to_order[m_coord] - 1;
                body.append_op(
                    "DETECTOR", {(k + 1) | TARGET_RECORD_BIT, (k + 1 + m) | TARGET_RECORD_BIT}, {m_coord.x, m_coord.y, 0});
            }
        } else {
            for (auto measure : chosen_basis_measure_coords) {
                auto k = (uint32_t)measurement_qubits.size() - measure_coord_to_order[measure] - 1;
                body.append_op(
                    "DETECTOR", {(k + 1) | TARGET_RECORD_BIT, (k + 1 + m) | TARGET_RECORD_BIT}, {measure.x, measure.y, 0});
            }
        }
    }

    // Build the end of the circuit, getting out of the cycle state and terminating.
    // In particular, the data measurements create detectors that have to be handled special.
    // Also, the tail is responsible for identifying the logical observable.
    Circuit tail;
    if (params.use_swap_lru) {
        // Two final sets of CNOTs to reswap data and ancilla qubits.
        std::array<std::vector<uint32_t>, 2> final_cnot_targets;
        size_t prev_cycle = (params.rounds % lru_cycles) - 1;
        for (auto measure : z_measure_coords) {
            uint32_t pm = p2q[measure];
            if (swap_order[pm][prev_cycle] != INVALID_SWAP && swap_order[pm][prev_cycle] != pm) {
                uint32_t pd = swap_order[pm][prev_cycle];
                final_cnot_targets[0].push_back(pm);
                final_cnot_targets[0].push_back(pd);
                final_cnot_targets[1].push_back(pd);
                final_cnot_targets[1].push_back(pm);
            }
        }

        // Ignore errors -- in practice, you would just measure
        // the swapped qubits instead of swapping them back.
        if (!final_cnot_targets[0].empty()) {
            tail.append_op("TICK", {});
#ifdef USE_SWAPS
            tail.append_op("SWAP", final_cnot_targets[0]);
#else
            tail.append_op("CNOT", final_cnot_targets[0]);
            tail.append_op("TICK", {});
            tail.append_op("CNOT", final_cnot_targets[1]);
#endif
        }
    }

    params.append_measure(tail, data_qubits, "ZX"[is_memory_x]);
    // Detectors.
    for (auto measure : chosen_basis_measure_coords) {
        std::vector<uint32_t> detectors;
        for (auto delta : z_order) {
            auto data = measure + delta;
            if (p2q.find(data) != p2q.end()) {
                detectors.push_back((data_qubits.size() - data_coord_to_order[data]) | TARGET_RECORD_BIT);
            }
        }
        detectors.push_back(
            (data_qubits.size() + measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT);
        std::sort(detectors.begin(), detectors.end());
        tail.append_op("DETECTOR", detectors, {measure.x, measure.y, 1});
    }
    // Logical observable.
    std::vector<uint32_t> obs_inc;
    for (auto q : chosen_basis_observable) {
        obs_inc.push_back((data_qubits.size() - data_coord_to_order[q]) | TARGET_RECORD_BIT);
    }
    std::sort(obs_inc.begin(), obs_inc.end());
    tail.append_op("OBSERVABLE_INCLUDE", obs_inc, 0);

    // Combine to form final circuit.
    Circuit main_body;
    if (params.use_swap_lru) {
        Circuit repeated_part, remaining_part;
        if (params.rounds - 1 >= lru_cycles) {
            repeated_part = circuit_bodies[1] + circuit_bodies[2] + circuit_bodies[3];
            if (params.swap_lru_with_no_swap) {
                repeated_part += circuit_bodies[4];
            }
            repeated_part += circuit_bodies[0];
            repeated_part *= (params.rounds - 1) / lru_cycles;

        }
        size_t cycles_left = (params.rounds - 1) % lru_cycles;
        size_t c = 1;
        while (cycles_left--) {
            remaining_part += circuit_bodies[c++];
        }
        main_body = repeated_part + remaining_part;
    } else {
        main_body = circuit_bodies[0] * (params.rounds - 1);
    }
    Circuit full_circuit = head + main_body + tail;

    // Produce a 2d layout.
    std::map<std::pair<uint32_t, uint32_t>, std::pair<std::string, uint32_t>> layout;
    float scale = x_order[0].x == 0.5 ? 2 : 1;
    for (auto q : data_coords) {
        layout[{(uint32_t)(q.x * scale), (uint32_t)(q.y * scale)}] = {"d", p2q[q]};
    }
    for (auto q : x_measure_coords) {
        layout[{(uint32_t)(q.x * scale), (uint32_t)(q.y * scale)}] = {"X", p2q[q]};
    }
    for (auto q : z_measure_coords) {
        layout[{(uint32_t)(q.x * scale), (uint32_t)(q.y * scale)}] = {"Z", p2q[q]};
    }
    for (auto q : chosen_basis_observable) {
        layout[{(uint32_t)(q.x * scale), (uint32_t)(q.y * scale)}].first = "L";
    }

    return {
        full_circuit,
        layout,
        "# Legend:\n"
        "#     d# = data qubit\n"
        "#     L# = data qubit with logical observable crossing\n"
        "#     X# = measurement qubit (X stabilizer)\n"
        "#     Z# = measurement qubit (Z stabilizer)\n"};
}

GeneratedCircuit _generate_rotated_surface_code_circuit(const CircuitGenParameters &params, bool is_memory_x) {
    uint32_t d = params.distance;

    // Place data qubits.
    std::set<surface_coord> data_coords;
    std::vector<surface_coord> x_observable;
    std::vector<surface_coord> z_observable;
    for (float x = 0.5; x <= d; x++) {
        for (float y = 0.5; y <= d; y++) {
            surface_coord q{x * 2, y * 2};
            data_coords.insert(q);
            if (y == 0.5) {
                z_observable.push_back(q);
            }
            if (x == 0.5) {
                x_observable.push_back(q);
            }
        }
    }

    // Place measurement qubits.
    std::set<surface_coord> x_measure_coords;
    std::set<surface_coord> z_measure_coords;
    for (size_t x = 0; x <= d; x++) {
        for (size_t y = 0; y <= d; y++) {
            surface_coord q{(float)x * 2, (float)y * 2};
            bool on_boundary_1 = x == 0 || x == d;
            bool on_boundary_2 = y == 0 || y == d;
            bool parity = x % 2 != y % 2;
            if (on_boundary_1 && parity) {
                continue;
            }
            if (on_boundary_2 && !parity) {
                continue;
            }
            if (parity) {
                x_measure_coords.insert(q);
            } else {
                z_measure_coords.insert(q);
            }
        }
    }

    // Define interaction orders so that hook errors run against the error grain instead of with it.
    std::vector<surface_coord> z_order{
        {1, 1},
        {1, -1},
        {-1, 1},
        {-1, -1},
    };
    std::vector<surface_coord> x_order{
        {1, 1},
        {-1, 1},
        {1, -1},
        {-1, -1},
    };

    // Delegate.
    return _finish_surface_code_circuit(
        [&](surface_coord q) {
            q = q - surface_coord{0, fmodf(q.x, 2)};
            return (uint32_t)(q.x + q.y * (d + 0.5));
        },
        data_coords,
        x_measure_coords,
        z_measure_coords,
        params,
        x_order,
        z_order,
        x_observable,
        z_observable,
        is_memory_x);
}

GeneratedCircuit _generate_unrotated_surface_code_circuit(const CircuitGenParameters &params, bool is_memory_x) {
    uint32_t d = params.distance;
    assert(params.rounds > 0);

    // Place qubits.
    std::set<surface_coord> data_coords;
    std::set<surface_coord> x_measure_coords;
    std::set<surface_coord> z_measure_coords;
    std::vector<surface_coord> x_observable;
    std::vector<surface_coord> z_observable;
    for (size_t x = 0; x < 2 * d - 1; x++) {
        for (size_t y = 0; y < 2 * d - 1; y++) {
            surface_coord q{(float)x, (float)y};
            bool parity = x % 2 != y % 2;
            if (parity) {
                if (x % 2 == 0) {
                    z_measure_coords.insert(q);
                } else {
                    x_measure_coords.insert(q);
                }
            } else {
                data_coords.insert(q);
                if (x == 0) {
                    x_observable.push_back(q);
                }
                if (y == 0) {
                    z_observable.push_back(q);
                }
            }
        }
    }

    // Define interaction order. Doesn't matter so much for unrotated.
    std::vector<surface_coord> order{
        {1, 0},
        {0, 1},
        {0, -1},
        {-1, 0},
    };

    // Delegate.
    return _finish_surface_code_circuit(
        [&](surface_coord q) {
            return (uint32_t)(q.x + q.y * (2 * d - 1));
        },
        data_coords,
        x_measure_coords,
        z_measure_coords,
        params,
        order,
        order,
        x_observable,
        z_observable,
        is_memory_x);
}

GeneratedCircuit stim::generate_surface_code_circuit(const CircuitGenParameters &params) {
    if (params.task == "rotated_memory_x") {
        return _generate_rotated_surface_code_circuit(params, true);
    } else if (params.task == "rotated_memory_z") {
        return _generate_rotated_surface_code_circuit(params, false);
    } else if (params.task == "unrotated_memory_x") {
        return _generate_unrotated_surface_code_circuit(params, true);
    } else if (params.task == "unrotated_memory_z") {
        return _generate_unrotated_surface_code_circuit(params, false);
    } else {
        throw std::invalid_argument(
            "Unrecognized task '" + params.task +
            "'. Known surface_code tasks:\n"
            "    'rotated_memory_x': Initialize logical |+> in rotated code, protect with parity measurements, measure "
            "logical X.\n"
            "    'rotated_memory_z': Initialize logical |0> in rotated code, protect with parity measurements, measure "
            "logical Z.\n"
            "    'unrotated_memory_x': Initialize logical |+> in unrotated code, protect with parity measurements, "
            "measure logical X.\n"
            "    'unrotated_memory_z': Initialize logical |0> in unrotated code, protect with parity measurements, "
            "measure logical Z.\n"
            "");
    }
}
