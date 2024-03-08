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
    decoding_graph(std::make_unique<DecodingGraph>(circuit, flips_per_error)),
    detectors(),
    flags(),
    flag_edges()
{}

void
MatchingBase::load_syndrome(
        stim::simd_bits_range_ref<SIMD_WIDTH> syndrome,
        int c1,
        int c2,
        bool recompute_flags) 
{
    std::vector<uint64_t> all_dets = get_nonzero_detectors(syndrome);
    detectors.clear();
    if (recompute_flags) {
        flags.clear();
        flag_edges.clear();
    }

#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;

    timer.clk_start();
#endif

    // Track the number of detectors of each color.
    for (auto it = all_dets.begin(); it != all_dets.end(); ) {
        uint64_t d = *it;
        if (circuit.flag_detectors.count(d)) {
            if (recompute_flags) {
                flags.push_back(d);
            }
            it = all_dets.erase(it);
        } else {
            if (c1 != COLOR_ANY && circuit.detector_color_map[d] != c1 && circuit.detector_color_map[d] != c2) {
                it = all_dets.erase(it);
            } else {
                it++;
            }
        }
    }
    detectors = std::move(all_dets);
    if (detectors.size() & 1) {
        detectors.push_back(get_color_boundary_index(COLOR_ANY));
    }
    // Activate flag edges in decoding_graph.
    if (recompute_flags) {
        decoding_graph->activate_flags(flags);
        flag_edges = decoding_graph->get_flag_edges();

        /*
        std::cout << "Flag edges:" << std::endl;
        for (sptr<hyperedge_t> e : flag_edges) {
            std::cout << "\tD[";
            for (size_t i = 0; i < e->get_order(); i++) {
                std::cout << " " << print_v(e->get<vertex_t>(i));
            }
            std::cout << " ], F[";
            for (uint64_t f : e->flags) {
                std::cout << " " << f;
            }
            std::cout << " ], frames:";
            for (uint64_t fr : e->frames) {
                std::cout << " " << fr;
            }
            std::cout << std::endl;
        }
        */
    }
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ MatchingBase ] Took " << t*1e-9 << "s to retrieve syndrome: D[";
    for (uint64_t d : detectors) std::cout << " " << d;
    std::cout << " ], F[";
    for (uint64_t f : flags) std::cout << " " << f;
    std::cout << " ]" << std::endl;
#endif
}

stim::simd_bits<SIMD_WIDTH>
MatchingBase::get_base_corr() {
#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;

    timer.clk_start();
#endif

    stim::simd_bits<SIMD_WIDTH> corr(circuit.count_observables());
    if (flags.empty()) {
#ifdef DECODER_PERF
        t = timer.clk_end();
        std::cout << "[ MatchingBase ] took " << t*1e-9 << "s to compute base correction" << std::endl;
#endif
        return corr;
    }
    // Otherwise, get a no-detector edge. If such an edge exists, then return the
    // corresponding correction.
    sptr<hyperedge_t> e = decoding_graph->get_best_nod_edge(detectors.size());
    if (e != nullptr) {
        for (uint64_t fr : e->frames) corr[fr] ^= 1;
    }
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ MatchingBase ] took " << t*1e-9 << "s to compute base correction" << std::endl;
#endif
    return corr;
}

Decoder::result_t
MatchingBase::ret_no_detectors() {
    return { 0.0, get_base_corr() };
}

std::vector<Decoder::assign_t>
MatchingBase::compute_matching(int c1, int c2, bool split_thru_boundary_match) {
    const size_t n = detectors.size();
    const size_t m = n*(n+1)/2;

    PerfectMatching pm(n, m);
    pm.options.verbose = false;
    // Add edges:
    std::map<uint64_t, uint64_t> boundary_pref_map;

#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;

    timer.clk_start();
#endif

    if (flags.size() && detectors.size() > 10) {
        decoding_graph->immediately_initialize_distances_for(c1, c2);
    }

    for (size_t i = 0; i < n; i++) {
        uint64_t di = detectors.at(i);
        for (size_t j = i+1; j < n; j++) {
            uint64_t dj = detectors.at(j);

            fp_t w;
            if (dj == get_color_boundary_index(COLOR_ANY) && c1 != COLOR_ANY) {
                // We need to identify the best boundary for di and use that.
                uint64_t b1 = get_color_boundary_index(c1),
                         b2 = get_color_boundary_index(c2);
                error_chain_t ec1 = decoding_graph->get_error_chain(di, b1, c1, c2),
                              ec2 = decoding_graph->get_error_chain(di, b2, c1, c2);
                boundary_pref_map[di] = ec1.weight < ec2.weight ? b1 : b2;
                w = std::min(ec1.weight, ec2.weight);
            } else {
                error_chain_t ec = decoding_graph->get_error_chain(di, dj, c1, c2);
                w = ec.weight;
            }
            uint32_t iw = static_cast<uint32_t>(1000 * w);
            pm.AddEdge(i, j, iw);
        }
    }
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ MatchingBase ] took " << t*1e-9 << "s to accumulate matching edges." << std::endl;
#endif

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
        error_chain_t ec = decoding_graph->get_error_chain(di, dj, c1, c2);
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

sptr<hyperedge_t>
MatchingBase::get_flag_edge_for(std::vector<sptr<vertex_t>> vlist) {
    for (sptr<hyperedge_t> e : flag_edges) {
        if (e->get_order() != vlist.size()) continue;

        bool all_match = true;
        for (size_t i = 0; i < e->get_order(); i++) {
            sptr<vertex_t> v = e->get<vertex_t>(i);
            if (std::find(vlist.begin(), vlist.end(), v) == vlist.end()) {
                all_match = false;
                break;
            }
        }
        if (all_match) return e;
    }
    return nullptr;
}

}   // qontra
