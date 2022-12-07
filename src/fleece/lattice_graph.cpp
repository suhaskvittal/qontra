/*
 *  author: Suhas Vittal
 *  date:   7 December 2022
 * */

#include "fleece/lattice_graph.h"

namespace qrc {
namespace fleece {

static LatticeGraph::Vertex NULL_VERTEX(-1, false, -1);

LatticeGraph::LatticeGraph()
    :vertex_list(),
    qubit_to_vertex(),
    detector_to_vertex(),
    adjacency_list()
{}

void
LatticeGraph::add_qubit(int32_t qubit, bool is_data, int32_t base_detector, int32_t meas_time) {
    if (qubit_to_vertex.count(qubit)) {
        // Update the data.
        qubit_to_vertex[qubit].is_data = is_data;
        qubit_to_vertex[qubit].base_detector = base_detector;
        if (meas_time >= 0) {
            qubit_to_vertex[qubit].measurement_times.push_back(meas_time);
        }
        if (base_detector >= 0) {
            detector_to_vertex[base_detector] = qubit_to_vertex[qubit];
        }
        for (Vertex& v : vertex_list) {
            if (v.qubit == qubit) {
                v.is_data = is_data;
                v.base_detector = base_detector;
                if (meas_time >= 0) {
                    v.measurement_times.push_back(meas_time);
                }
            }
        }
        return;
    }

    Vertex v(qubit, is_data, base_detector);
    vertex_list.push_back(v);
    qubit_to_vertex[qubit] = v;
    if (base_detector >= 0) {
        detector_to_vertex[base_detector] = v;
    }
    if (meas_time >= 0) {
        v.measurement_times.push_back(meas_time);
    }
    adjacency_list[v] = std::vector<Vertex>();
}

void
LatticeGraph::add_coupling(int32_t q1, int32_t q2) {
    Vertex v1 = qubit_to_vertex[q1];
    Vertex v2 = qubit_to_vertex[q2];
    add_coupling(v1, v2);
}

void
LatticeGraph::add_coupling(const Vertex& v1, const Vertex& v2) {
    // Check that we aren't causing a duplicate.
    for (const Vertex& w : adjacency_list[v1]) {
        if (w == v2) {
            return;
        }
    } 
    
    adjacency_list[v1].push_back(v2);
    adjacency_list[v2].push_back(v1);
}

LatticeGraph::Vertex
LatticeGraph::get_vertex_by_qubit(int32_t q) {
    if (qubit_to_vertex.count(q)) {
        return qubit_to_vertex[q];
    } else {
        return NULL_VERTEX;
    }
}

LatticeGraph::Vertex
LatticeGraph::get_vertex_by_detector(int32_t d) {
    if (detector_to_vertex.count(d)) {
        return detector_to_vertex[d];
    } else {
        return NULL_VERTEX;
    }
}

std::vector<LatticeGraph::Vertex>
LatticeGraph::vertices() {
    return vertex_list;
}

std::vector<LatticeGraph::Vertex>
LatticeGraph::get_adjacency_list(int32_t q) {
    Vertex v = qubit_to_vertex[q];
    return get_adjacency_list(v);
}

std::vector<LatticeGraph::Vertex>
LatticeGraph::get_adjacency_list(const Vertex& v) {
    return adjacency_list[v];
}

LatticeGraph
to_lattice_graph(const stim::Circuit& circuit) {
    std::deque<int32_t> measurement_order;
    std::set<int32_t> already_measured_once;

    uint32_t detector_counter = 0;

    LatticeGraph graph;

    for (const stim::Operation& op : circuit.operations) {
        std::string opname(op.gate->name);  
        if (opname == "QUBIT_COORDS") {
            // This is a declaration of a qubit. Create a vertex.
            int32_t qubit = (int32_t)op.target_data.targets[0].data;
            graph.add_qubit(qubit, true, -1);
        } else if (opname == "CX" || opname == "ZCX") {
            // This is a CNOT, indicating a coupling.
            const auto& targets = op.target_data.targets;
            for (uint32_t i = 0; i < targets.size(); i += 2) {
                int32_t q1 = (int32_t)targets[i].data;
                int32_t q2 = (int32_t)targets[i+1].data;
                graph.add_coupling(q1, q2);
            }
        } else if (opname == "M" || opname == "MX" || opname == "MZ" || opname == "MR") {
            const auto& targets = op.target_data.targets;
            for (auto target : targets) {
                graph.add_qubit((int32_t)target.data, false, -1, measurement_order.size());
                measurement_order.push_front((int32_t)target.data);
            } 
        } else if (opname == "DETECTOR") {
            const auto& targets = op.target_data.targets;
            for (auto target : targets) {
                int32_t q = measurement_order[(int32_t)target.data];
                if (!already_measured_once.count(q)) {
                    graph.add_qubit(q, false, detector_counter++);
                    already_measured_once.insert(q);
                }
            }
        } else if (opname == "TAILSTART") {
            break;
        }
    }

    return graph;
}

}  // fleece
}  // qrc
