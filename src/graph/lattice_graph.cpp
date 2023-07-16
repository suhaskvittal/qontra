/*
 *  author: Suhas Vittal
 *  date:   16 July 2023
 * */

#include "graph/lattice_graph.h"

namespace qontra {
namespace graph {

using namespace lattice;

bool
LatticeGraph::update_state() {
    if (!__LatticeGraphParent::update_state())  return false;
    build_distance_matrix();
    return true;
}

void
LatticeGraph::build_distance_matrix() {
    distance::callback_t<vertex_t, matrix_entry_t> d_cb =
        [&] (vertex_t* src,
                vertex_t* dst,
                const std::map<vertex_t*, fp_t>& dist,
                const std::map<vertex_t*, vertex_t*>& pred)
        {
            matrix_entry_t entry;

            auto curr = dst;
            while (curr != src) {
                if (curr->qubit_type == vertex_t::DATA) {
                    entry.error_chain.push_back(curr);
                }
                entry.physical_path.push_back(curr);
                curr = pred.at(curr);
            }

            std::reverse(entry.error_chain.begin(), entry.error_chain.end());
            std::reverse(entry.physical_path.begin(), entry.physical_path.end());

            return entry;
        };
    distance_matrix = distance::create_distance_matrix(
                            this, unit_ewf_t<vertex_t>(), d_cb);
}

namespace io {

void
update_lattice_graph(LatticeGraph& lattice_graph, std::string line) {
    size_t ssi = line.find(",");
    
    std::string line_type = line.substr(0, ssi);
    // Now, get the data from the line.
    size_t pssi = ssi+1;
    
    uint32_t edge_src_id;
    bool waiting_for_edge_dst = false;
    while (true) {
        if ((ssi = line.find(",", pssi)) == std::string::npos) {
            ssi = line.size();  // We have hit the end of the line.
        }
        uint32_t n = std::stoi(line.substr(pssi, ssi-pssi));
        // Make sure the graph contains the vertex before doing anything.
        if (!lattice_graph.contains(n)) {
            vertex_t* v = new vertex_t;
            v->id = n;
            lattice_graph.add_vertex(v);
        }
        vertex_t* v = lattice_graph.get_vertex(n);

        if (line_type == "D") {
            v->qubit_type = vertex_t::DATA;
        } else if (line_type == "PX") {
            v->qubit_type = vertex_t::XPARITY;
            lattice_graph.add_to_meas_order(v);
        } else if (line_type == "PZ") {
            v->qubit_type = vertex_t::ZPARITY;
            lattice_graph.add_to_meas_order(v);
        } else if (line_type == "E") {
            if (waiting_for_edge_dst) {
                vertex_t* w = lattice_graph.get_vertex(edge_src_id);
                if (!lattice_graph.contains(v, w)) {
                    auto e = new edge_t;
                    e->src = w;
                    e->dst = v;
                    e->is_undirected = true;
                    lattice_graph.add_edge(e);
                }
                waiting_for_edge_dst = false;
            } else {
                edge_src_id = n;
                waiting_for_edge_dst = true;
            }
        }

        if (ssi == line.size()) break;
        pssi = ssi+1;
    }
}

}   // io

}   // graph
}   // qontra
