/*
 *  author: Suhas Vittal
 *  date:   16 February 2024
 * */

#define MEMORY_DEBUG

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
{
    if (decoding_graph->number_of_colors == 0) {
        decoding_graph->immediately_initialize_distances_for(COLOR_ANY, COLOR_ANY);
    } else {
        for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
            for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
                decoding_graph->immediately_initialize_distances_for(c1, c2);
            }
        }
    }
}

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
        decoding_graph->activate_detectors(detectors, flags);
        flag_edges = decoding_graph->get_flag_edges();

#ifdef MEMORY_DEBUG
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
#endif
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

std::vector<assign_t>
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
    std::vector<assign_t> assign_arr;
    for (size_t i = 0; i < n; i++) {
        size_t j = pm.GetMatch(i);
        uint64_t di = detectors.at(i),
                 dj = detectors.at(j);
        if (di > dj) continue;
        // Replace the boundary if necessary.
        if (dj == get_color_boundary_index(COLOR_ANY) && c1 != COLOR_ANY) {
            dj = boundary_pref_map[di];
        }
        sptr<vertex_t> v = decoding_graph->get_vertex(di),
                        w = decoding_graph->get_vertex(dj);
        expand_error_chain(assign_arr, v, w, c1, c2, split_thru_boundary_match);
    }
    return assign_arr;
}

void
MatchingBase::expand_error_chain(
        std::vector<assign_t>& assign_arr,
        sptr<vertex_t> src,
        sptr<vertex_t> dst, 
        int c1, int c2,
        bool split_through_boundary_match) 
{
    error_chain_t ec = decoding_graph->get_error_chain(src, dst, c1, c2);
    
    bool all_entries_are_boundaries = ec.path.front()->is_boundary_vertex;
    assign_t curr_assign;
    curr_assign.v = ec.path.front();
    curr_assign.c1 = c1;
    curr_assign.c2 = c2;
    curr_assign.path = {ec.path.front()};

    for (size_t i = 1; i < ec.path.size(); i++) {
        sptr<vertex_t> v = ec.path[i-1],
                        w = ec.path[i];
        if (decoding_graph->share_hyperedge({v, w})) {
            curr_assign.w = w;
            curr_assign.path.push_back(w);
            all_entries_are_boundaries &= w->is_boundary_vertex;
        } else {
            // This may be a flag edge.
            sptr<hyperedge_t> e = get_flag_edge_for({v, w});
            if (e == nullptr) {
#ifdef MEMORY_DEBUG
                std::cout << "\tForced to expand " << print_v(v) << ", " << print_v(w) << std::endl;
#endif
                // Apparently not. Fill in the gap.
                error_chain_t _ec = decoding_graph->get_error_chain(v, w, c1, c2, true);
                if (_ec.path.front() == w) std::reverse(_ec.path.begin(), _ec.path.end());
                for (size_t j = 1; j < _ec.path.size(); j++) {
                    sptr<vertex_t> x = _ec.path[j];
                    curr_assign.w = x;
                    curr_assign.path.push_back(x);
                    all_entries_are_boundaries &= x->is_boundary_vertex;
                    if (x->is_boundary_vertex && split_through_boundary_match) {
                        if (!all_entries_are_boundaries) {
                            assign_arr.push_back(curr_assign);
#ifdef MEMORY_DEBUG
                            std::cout << "\t[B2] Found assignment " << print_v(curr_assign.v) << " <---> "
                                << print_v(curr_assign.w) << " with path:";
                            for (sptr<vertex_t> x : curr_assign.path) std::cout << " " << print_v(x);
                            std::cout << std::endl;
#endif
                        } 
                        curr_assign.v = x;
                        curr_assign.w = nullptr;
                        curr_assign.path = {x};
                        curr_assign.flag_edges.clear();
                        all_entries_are_boundaries = true;
                    }
                }
            } else {
                curr_assign.flag_edges.push_back(e);
                // Break the path here.
                if (curr_assign.w == nullptr) {
                    curr_assign.v = nullptr;
                }
                assign_arr.push_back(curr_assign);
#ifdef MEMORY_DEBUG
                std::cout << "\t[F] Found assignment " << print_v(curr_assign.v) << " <---> "
                    << print_v(curr_assign.w) << " with path:";
                for (sptr<vertex_t> x : curr_assign.path) std::cout << " " << print_v(x);
                std::cout << std::endl;
#endif
                curr_assign.v = w;
                curr_assign.w = nullptr;
                curr_assign.path = {w};
                curr_assign.flag_edges.clear();
                all_entries_are_boundaries = w->is_boundary_vertex;
                continue;
            }
        }
        // Now handle any boundaries if necessary.
        if (!split_through_boundary_match || !w->is_boundary_vertex) continue;
        if (!all_entries_are_boundaries) {
            assign_arr.push_back(curr_assign);
#ifdef MEMORY_DEBUG
            std::cout << "\t[B] Found assignment " << print_v(curr_assign.v) << " <---> "
                << print_v(curr_assign.w) << " with path:";
            for (sptr<vertex_t> x : curr_assign.path) std::cout << " " << print_v(x);
            std::cout << std::endl;
#endif
        }
        curr_assign.v = w;
        curr_assign.w = nullptr;
        curr_assign.path = {w};
        curr_assign.flag_edges.clear();
        all_entries_are_boundaries = true;
    }
    // Handle the remaining assignment.
    if (!all_entries_are_boundaries && curr_assign.w != nullptr) {
#ifdef MEMORY_DEBUG
        std::cout << "\t[T] Found assignment " << print_v(curr_assign.v) << " <---> "
            << print_v(curr_assign.w) << " with path:";
        for (sptr<vertex_t> x : curr_assign.path) std::cout << " " << print_v(x);
        std::cout << std::endl;
#endif
        assign_arr.push_back(curr_assign);
    }
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
