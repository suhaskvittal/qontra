/*
 *  author: Suhas Vittal
 *  date:   7 December 2022
 * */

#include "fleece/lattice_graph.h"

namespace qrc {
namespace fleece {

LatticeGraph::LatticeGraph()
    :vertex_list(),
    qubit_to_vertex(),
    detector_to_vertex(),
    adjacency_matrix()
{}

LatticeGraph::~LatticeGraph() {
}

bool
LatticeGraph::add_qubit(int32_t qubit, uint8_t is_data,
        int32_t detector, int32_t meas_time, bool force_record) 
{
    if (qubit_to_vertex.count(qubit)) {
        // Update the data.
        if (!(is_data & 0b10)) {
            qubit_to_vertex[qubit]->is_data = is_data;
        }
        if (detector >= 0 || force_record) {
            qubit_to_vertex[qubit]->detectors.push_back(detector);
            detector_to_vertex[detector] = qubit_to_vertex[qubit];
        }
        if (meas_time >= 0 || force_record) {
            qubit_to_vertex[qubit]->measurement_times.push_back(meas_time);
        }
        return false;
    }

    Vertex * v = new Vertex(qubit, is_data);
    vertex_list.push_back(v);
    qubit_to_vertex[qubit] = v;
    if (detector >= 0 || force_record) {
        v->detectors.push_back(detector);
        detector_to_vertex[detector] = v;
    }
    if (meas_time >= 0 || force_record) {
        v->measurement_times.push_back(meas_time);
    }
    adjacency_matrix[v] = std::vector<Vertex*>();
    return true;
}

bool
LatticeGraph::add_coupling(int32_t q1, int32_t q2) {
    return add_coupling(qubit_to_vertex[q1], qubit_to_vertex[q2]);
}

bool
LatticeGraph::add_coupling(Vertex * v1, Vertex * v2) {
    // Check that we aren't causing a duplicate.
    for (Vertex * w : adjacency_matrix[v1]) {
        if (w == v2) {
            return false;
        }
    } 
    
    adjacency_matrix[v1].push_back(v2);
    adjacency_matrix[v2].push_back(v1);
    return true;
}

LatticeGraph::Vertex*
LatticeGraph::get_vertex_by_qubit(int32_t q) {
    if (qubit_to_vertex.count(q)) {
        return qubit_to_vertex[q];
    } else {
        return nullptr;
    }
}

LatticeGraph::Vertex*
LatticeGraph::get_vertex_by_detector(int32_t d) {
    if (detector_to_vertex.count(d)) {
        return detector_to_vertex[d];
    } else {
        return nullptr;
    }
}

std::vector<LatticeGraph::Vertex*>
LatticeGraph::vertices() {
    return vertex_list;
}

std::vector<LatticeGraph::Vertex*>
LatticeGraph::adjacency_list(int32_t q) {
    return adjacency_list(qubit_to_vertex[q]);
}

std::vector<LatticeGraph::Vertex*>
LatticeGraph::adjacency_list(Vertex * v) {
    return adjacency_matrix[v];
}

LatticeGraph
to_lattice_graph(const stim::Circuit& circuit) {
    std::deque<int32_t> measurement_order;
    std::map<int32_t, uint8_t> neighbor_count;

    bool first_round = true;
    uint32_t detector_counter = 0;
    uint32_t measurement_time = 0;

    LatticeGraph graph;
    stim::Circuit flat_circ = circuit.flattened();

    for (const stim::Operation& op : flat_circ.operations) {
        std::string opname(op.gate->name);  
        if (opname == "QUBIT_COORDS") {
            // This is a declaration of a qubit. Create a vertex.
            int32_t qubit = (int32_t)op.target_data.targets[0].data;
            graph.add_qubit(qubit, true, -1);
        } else if (opname == "CX" || opname == "ZCX" && first_round) {
            // This is a CNOT, indicating a coupling.
            const auto& targets = op.target_data.targets;
            for (uint32_t i = 0; i < targets.size(); i += 2) {
                int32_t q1 = (int32_t)targets[i].data;
                int32_t q2 = (int32_t)targets[i+1].data;
                if (graph.add_coupling(q1, q2)) {
                    if (!neighbor_count.count(q1)) {
                        neighbor_count[q1] = 0;
                    }
                    if (!neighbor_count.count(q2)) {
                        neighbor_count[q2] = 0;
                    }
                    neighbor_count[q1]++;
                    neighbor_count[q2]++;
                }
            }
        } else if (opname == "M" || opname == "MX" 
                || opname == "MZ" || opname == "MR") 
        {
            const auto& targets = op.target_data.targets;
            for (auto target : targets) {
                uint8_t is_data = (!first_round << 1);
                graph.add_qubit(
                        (int32_t)target.data, is_data, -1, measurement_time++);
                measurement_order.push_front((int32_t)target.data);
            } 
        } else if (opname == "DETECTOR") {
            const auto& targets = op.target_data.targets;
            for (auto target : targets) {
                int32_t index = (int32_t)(target.data ^ stim::TARGET_RECORD_BIT) - 1;
                if (index >= measurement_order.size()) {
                    continue;
                }
                int32_t q = measurement_order[index];
                uint8_t is_data = (!first_round << 1);
                graph.add_qubit(q, is_data, detector_counter++);
            }
        } else if (opname == "SIMHALT") {
            first_round = false;
            measurement_order.clear(); 
        } else if (opname == "TAILSTART") {
            break;
        }
    }
    // Now, identify any data qubits connected to less than 4 parity qubits.
    // Connect these to the boundary.
    graph.add_qubit(LATTICE_BOUNDARY, false, LATTICE_BOUNDARY_DETECTOR, -1, true);
    for (auto pair : neighbor_count) {
        if (pair.second < 4) {
            graph.add_coupling(LATTICE_BOUNDARY, pair.first);
        }
    }

    return graph;
}

graph::PathTable<LatticeGraph::Vertex>
compute_path_table(LatticeGraph& graph) {
    typedef LatticeGraph G;
    typedef LatticeGraph::Vertex V;

    graph::ewf_t<G, V> w = [] (G& g, V * v1, V * v2)
    {
        return 1;
    };

    return graph::compute_path_table(graph, w);
}

}  // fleece
}  // qrc
