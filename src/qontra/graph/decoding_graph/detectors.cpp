/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"

namespace qontra {
namespace graph {

using namespace decoding;

void
DecodingGraph::activate_detectors(const std::vector<uint64_t>& all_detectors) {
    deactivate_detectors();
    for (uint64_t d : all_detectors) {
        if (flag_detectors.count(d)) {
            active_flags.insert(d);
        } else {
            active_detectors.insert(d);
        }
    }
    flags_are_active = reweigh_for_detectors || !active_flags.empty();
    renorm_factor = compute_renorm_factor();
}

void
DecodingGraph::activate_detectors(
        const std::vector<uint64_t>& nonflags, const std::vector<uint64_t>& flags) 
{
    deactivate_detectors();
    active_detectors = std::unordered_set<uint64_t>(nonflags.begin(), nonflags.end());
    active_flags = std::unordered_set<uint64_t>(flags.begin(), flags.end());
    flags_are_active = reweigh_for_detectors || !active_flags.empty();
    renorm_factor = compute_renorm_factor();
}

void
DecodingGraph::deactivate_detectors() {
    flagged_distance_matrix_map.clear();
    flagged_dijkstra_graph_map.clear();
    active_flags.clear();
    active_detectors.clear();
    flags_are_active = false;
}

}   // graph
}   // qontra
