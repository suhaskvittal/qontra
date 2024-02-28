/*
 *  author: Suhas Vittal
 *  date:   16 February 2024
 * */

#include "qontra/decoder/matching_base.h"

#include <PerfectMatching.h>

namespace qontra {

using namespace graph;
using namespace decoding;

MatchingBase::MatchingBase(const DetailedStimCircuit& circuit, int flips_per_error)
    :Decoder(circuit),
    base_corr(circuit.count_observables()),
    decoding_graph(std::make_unique<DecodingGraph>(circuit, flips_per_error)),
    detectors(),
    flags()
{}

void
MatchingBase::load_syndrome(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome, int c1, int c2) {
    base_corr.clear();

    std::vector<uint64_t> all_dets = get_nonzero_detectors(syndrome);
    detectors.clear();
    flags.clear();
    // Track the number of detectors of each color.
    for (uint64_t d : all_dets) {
        if (circuit.flag_detectors.count(d)) {
            flags.push_back(d);
        } else {
            if (c1 != COLOR_ANY && circuit.detector_color_map[d] != c1 && circuit.detector_color_map[d] != c2) {
                continue;
            }
            detectors.push_back(d);
        }
    }
    // Activate flag edges in decoding_graph.
    decoding_graph->activate_flags(flags);
    // Check for any flag singletons. If any exist, remove any detectors associated with the singleton.
    if (flags.size()) {
        for (sptr<hyperedge_t> e : decoding_graph->get_flag_singletons()) {
            uint64_t d = e->get<vertex_t>(0)->id;
            auto it = std::find(detectors.begin(), detectors.end(), d);
            if (it != detectors.end()) {
                detectors.erase(it);
                for (uint64_t fr : e->frames) base_corr[fr] ^= 1;
            }
        }
    }
    if (detectors.size() & 1) {
        detectors.push_back(get_color_boundary_index(COLOR_ANY));
    }
}

std::vector<Decoder::assign_t>
MatchingBase::compute_matching(int c1, int c2, bool split_thru_boundary_match) {
    const size_t n = detectors.size();
    const size_t m = n*(n+1)/2;

    PerfectMatching pm(n, m);
    pm.options.verbose = false;
    // Add edges:
    std::map<uint64_t, uint64_t> boundary_pref_map;

    for (size_t i = 0; i < n; i++) {
        uint64_t di = detectors.at(i);
        for (size_t j = i+1; j < n; j++) {
            uint64_t dj = detectors.at(j);

            fp_t w;
            if (dj == get_color_boundary_index(COLOR_ANY) && c1 != COLOR_ANY) {
                // We need to identify the best boundary for di and use that.
                uint64_t b1 = get_color_boundary_index(c1),
                         b2 = get_color_boundary_index(c2);
                error_chain_t ec1 = decoding_graph->get(c1, c2, di, b1),
                              ec2 = decoding_graph->get(c1, c2, di, b2);
                boundary_pref_map[di] = ec1.weight < ec2.weight ? b1 : b2;
                w = std::min(ec1.weight, ec2.weight);
            } else {
                error_chain_t ec = decoding_graph->get(c1, c2, di, dj);
                w = ec.weight;
                // Quantize the weight.
            }
            uint32_t iw = w > 1000.0 ? 1'000'000 : static_cast<uint32_t>(1000 * w);
            pm.AddEdge(i, j, iw);
        }
    }
    pm.Solve();
    // Return assignments. Deactivate flags to avoid using flag edges.
    std::vector<Decoder::assign_t> assign_arr;
    for (size_t i = 0; i < n; i++) {
        size_t j = pm.GetMatch(i);
        uint64_t di = detectors.at(i),
                 dj = detectors.at(j);
        if (di > dj) continue;
        // Replace the boundary if necessary.
        if (dj == get_color_boundary_index(COLOR_ANY) && c1 != COLOR_ANY) {
            dj = boundary_pref_map[di];
        }
        error_chain_t ec = decoding_graph->get(c1, c2, di, dj);
        if (split_thru_boundary_match && ec.runs_through_boundary) {
            // Partition path between di and dj with boundaries.
            bool all_entries_are_boundaries = true;
            std::vector<uint64_t> part;
            for (size_t i = 0; i < ec.path.size(); i++) {
                auto v = ec.path.at(i);
                uint64_t d = v->id;
                part.push_back(d);
                all_entries_are_boundaries &= v->is_boundary_vertex;
                if (v->is_boundary_vertex) {
                    if (!all_entries_are_boundaries) {
                        // Add the endpoints as an assignment.
                        uint64_t id1 = part[0],
                                    id2 = part.back();
                        if (id1 > id2) std::swap(id1, id2);
                        assign_arr.emplace_back(id1, id2);
                    }
                    part = { d };
                    all_entries_are_boundaries = true;
                }
            }
            // Add the remaining part as an assignment.
            if (!all_entries_are_boundaries) {
                uint64_t id1 = part[0],
                            id2 = part.back();
                if (id1 > id2) std::swap(id1, id2);
                assign_arr.emplace_back(id1, id2);
            }
        } else {
            assign_arr.emplace_back(di, dj);
        }
    }
    return assign_arr;
}

}   // qontra
