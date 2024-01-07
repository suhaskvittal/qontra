/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include <algorithm>

namespace qontra {
namespace graph {

inline LatticeGraph::matrix_entry_t
LatticeGraph::get_path_data(sptr<lattice::vertex_t> v1, sptr<lattice::vertex_t> v2) {
    update_state();
    return distance_matrix[v1][v2];
}

inline bool
LatticeGraph::update_state() {
    if (!__LatticeGraphParent::update_state())  return false;
    build_distance_matrix();
    return true;
}

inline void
LatticeGraph::build_distance_matrix() {
    using namespace lattice;
    distance_matrix = create_distance_matrix<vertex_t, edge_t, matrix_entry_t>(this, 
                        // Weight function
                        [] (sptr<edge_t> e) {
                            return 1.0;
                        }, 
                        // Dijkstra callback:
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
                        });
}

}   // graph
}   // qontra
