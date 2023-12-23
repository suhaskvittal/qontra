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
        [&] (sptr<vertex_t> src,
                sptr<vertex_t> dst,
                const std::map<sptr<vertex_t>, fp_t>& dist,
                const std::map<sptr<vertex_t>, sptr<vertex_t>>& pred)
        {
            matrix_entry_t entry;

            auto curr = dst;
            while (curr != src) {
                if (curr->qubit_type == vertex_t::type::data) {
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

}   // graph
}   // qontra
