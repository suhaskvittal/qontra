/*
 *  author: Suhas Vittal
 *  date:   16 February 2024
 * */

#include "qontra/decoder/matching_base.h"

#include <PerfectMatching.h>

namespace qontra {

using namespace graph;

void
MatchingBase::load_syndrome(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome, int c1, int c2) {
    std::vector<uint64_t> all_dets = get_nonzero_detectors(syndrome);
    detectors.clear();
    flags.clear();
    // Track the number of detectors of each color.
    size_t n_of_c1 = 0,
           n_of_c2 = 0;
    for (uint64_t d : all_dets) {
        if (circuit.flag_detectors.count(d)) {
            flags.push_back(d);
        } else {
            if (circuit.detector_color_map[d] == c1) {
                n_of_c1++;
            } else if (circuit_detector_color_map[d] == c2) {
                n_of_c2++;
            } else {
                continue;
            }
            detectors.push_back(d);
        }
    }
    if (detectors.size() & 1) {
        if (c1 == COLOR_ANY) {
            detectors.push_back(get_color_boundary_index(COLOR_ANY));
        } else if (n_of_c1 & 1) {
            // Push back a boundary of color c2.
            detectors.push_back(get_color_boundary_index(c2));
        } else {
            detectors.push_back(get_color_boundary_index(c1));
        }
    }
    // Activate flag edges in decoding_graph.
    decoding_graph->activate_flags(flags);
}

std::vector<Decoder::assign_t>
MatchingBase::compute_matching(int c1, int c2, bool split_thru_boundary_match) {
    const size_t n = detectors.size();
    const size_t m = n*(n+1)/2;

    PerfectMatching pm(n, m);
    pm.options.verbose = false;
    // Add edges:
    for (size_t i = 0; i < n; i++) {
        uint64_t di = detectors.at(i);
        for (size_t j = i+1; j < n; j++) {
            uint64_t dj = detectors.at(j);
            error_chain_t ec = decoding_graph->get(c1, c2, di, dj);
            // Quantize the weight.
            uint32_t iw = static_cast<uint32_t>(1024 * ec.weight);
            pm.AddEdge(i, j, iw);
        }
    }
    pm.Solve();
    // Return assignments.
    std::vector<Decoder::assign_t> assign_arr;
    for (size_t i = 0; i < n; i++) {
        size_t j = pm.GetMatch(i);
        if (i > j) continue;
        uint64_t di = detectors.at(i),
                 dj = detectors.at(j);
        error_chain_t ec = decoding_graph->get(c1, c2, di, dj);
        if (split_thru_boundary_match && ec.runs_through_boundary) {
            // Partition path between di and dj with boundaries.
            bool all_entries_are_boundaries = true;
            std::vector<uint64_t> part;
            for (size_t i = 0; i < ec.path.size(); i++) {
                auto v = ec.path.at(i);
                uint64_t d = v->id;
                part.push_back(d);
                all_entries_are_boundaries |= !(v->is_boundary);
                if (v->is_boundary) {
                    if (!all_entries_are_boundaries) {
                        // Add the endpoints as an assignment.
                        assign_arr.emplace_back(part[0], part.back());
                    }
                    part.clear();
                    all_entries_are_boundaries = true;
                }
            }
            // Add the remaining part as an assignment.
            if (!all_entries_are_boundaries) {
                assign_arr.emplace_back(part[0], part.back());
            }
        } else {
            assign_arr.emplace_back(di, dj);
        }
    }
    return assign_arr;
}

}   // qontra
