/*
 *  author: Suhas Vittal
 *  date:   27 May 2024
 * */

#include "qontra/decoder/concat_mwpm.h"

#include <vtils/utility.h>

#include <PerfectMatching.h>

#include <math.h>

namespace qontra {

using namespace graph;
using namespace decoding;

ConcatMWPMDecoder::ConcatMWPMDecoder(const DetailedStimCircuit& circuit)
    :MatchingBase(circuit, 3, false),
    rgb_only_graphs(),
    edge_vertex_map()
{
    for (int c = 0; c < 3; c++) {
        rgb_only_graphs[c] = std::make_unique<DecodingGraph>(
                                decoding_graph->make_rgb_only_lattice(c, edge_vertex_map));
        rgb_only_graphs[c]->init_distances_for(COLOR_ANY, COLOR_ANY);
    }
}

Decoder::result_t
ConcatMWPMDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();

    load_syndrome(syndrome);
    if (detectors.empty()) return ret_no_detectors();

#ifdef MEMORY_DEBUG
    std::cout << "Detectors:";
    for (uint64_t d : detectors) std::cout << " " << d;
    std::cout << ", Flags:";
    for (uint64_t f : flags) std::cout << " " << f;
    std::cout << std::endl;
#endif

    std::array<std::vector<sptr<vertex_t>>, 3> detectors_by_color;
    for (uint64_t d : detectors) {
        if (d == get_color_boundary_index(COLOR_ANY)) continue;
        sptr<vertex_t> v = decoding_graph->get_vertex(d);
        detectors_by_color[v->color].push_back(v);
    }
    
    fp_t best_log_pr = std::numeric_limits<fp_t>::lowest();
    stim::simd_bits<SIMD_WIDTH> best_corr(n_obs);
    for (int c = 0; c < 3; c++) {
#ifdef MEMORY_DEBUG
        std::cout << "COLOR = " << c << std::endl;
#endif
        int cr1 = (c+1)%3, cr2 = (c+2)%3;
        load_syndrome(syndrome, cr1, cr2, false);
        std::vector<assign_t> matchings = compute_matching(cr1,cr2);
        // Identify vertex-vertex and edge-vertex nodes from the matching.
        std::set<sptr<vertex_t>> _active_vertices;
        for (const assign_t& m : matchings) {
            for (size_t i = 1; i < m.path.size(); i++) {
                sptr<vertex_t> v = m.path.at(i-1),
                                w = m.path.at(i);
                if (!decoding_graph->share_hyperedge({v,w})) {
                    // This is a flag edge. Get path between v and w.
                    error_chain_t ec = decoding_graph->get_error_chain(v, w, cr1, cr2, true);
                    if (ec.path.front() != v) std::reverse(ec.path.begin(), ec.path.end());
                    for (size_t j = 1; j < ec.path.size(); j++) {
                        sptr<vertex_t> x = ec.path.at(j-1),
                                        y = ec.path.at(j);
                        if (x->is_boundary_vertex && y->is_boundary_vertex) continue;
                        auto xy = edge_vertex_map.at(make_ev_pair(x,y));
                        _active_vertices.insert(xy);
#ifdef MEMORY_DEBUG
                        std::cout << "\tMatching " << print_v(x)
                            << " and " << print_v(y) << " --> "
                            << print_v(xy) << std::endl;
#endif
                    }
                } else {
                    if (v->is_boundary_vertex && w->is_boundary_vertex) continue;
                    auto vw = edge_vertex_map.at(make_ev_pair(v,w));
                    _active_vertices.insert(vw);
#ifdef MEMORY_DEBUG
                    std::cout << "\tMatching " << print_v(v) << " and " << print_v(w)
                        << " --> " << print_v(vw) << std::endl;
#endif
                }
            }
        }
        if (_active_vertices.empty()) continue;
        std::vector<sptr<vertex_t>> 
            active_vertices(_active_vertices.begin(), _active_vertices.end());
        vtils::push_back_range(active_vertices, detectors_by_color[c]);
        if (active_vertices.size() & 1) {
            auto vb = decoding_graph->get_boundary_vertex(c);
            // Check if the boundary is a false boundary (occur in even distance codes)
            if (decoding_graph->get_degree(vb) <= 1) continue;
#ifdef MEMORY_DEBUG
            std::cout << "\tAdding boundary vertex." << std::endl;
#endif
            active_vertices.push_back(vb);
        }
        // Add red vertices as well.
        stim::simd_bits<SIMD_WIDTH> corr(n_obs);
        fp_t log_pr = rgb_compute_matching(active_vertices, c, corr);
        if (log_pr > best_log_pr) {
            best_log_pr = log_pr;
            best_corr = corr;
        }
    }
    return {0.0, best_corr};
}


fp_t
ConcatMWPMDecoder::rgb_compute_matching(
        const std::vector<sptr<vertex_t>>& active,
        int color,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr)
{
    uptr<DecodingGraph>& gr = rgb_only_graphs[color];

    std::vector<uint64_t> local_dets;
    for (auto v : active) {
        if (!v->is_boundary_vertex) {
            local_dets.push_back(v->id);
        }
    }
    gr->activate_detectors(local_dets, flags);
    // Now, create the matching graph.

    const size_t n = active.size();
    const size_t m = n*(n+1)/2;
    PerfectMatching pm(n,m);
    pm.options.verbose = false;

    for (size_t i = 0; i < n; i++) {
        auto v = active.at(i);
        for (size_t j = i+1; j < n; j++) {
            auto w = active.at(j);
            error_chain_t ec = gr->get_error_chain(v, w);
            uint32_t iw = static_cast<uint32_t>(1000 * ec.weight);
            pm.AddEdge(i, j, iw);
        }
    }
    pm.Solve();
    // Now compute the correction and its probability.
    fp_t log_pr = 0.0;
    for (size_t i = 0; i < n; i++) {
        size_t j = pm.GetMatch(i);
        if (i > j) continue;
        auto v = active.at(i),
             w = active.at(j);
        error_chain_t ec = gr->get_error_chain(v, w);
        // At this stage, flag edges do not matter.
#ifdef MEMORY_DEBUG
        std::cout << "\tMatched " << print_v(v) << " <---> " << print_v(w) 
            << ", probability = " << ec.probability << ", path size = "
            << ec.path.size() << std::endl;
#endif
        for (size_t k = 1; k < ec.path.size(); k++) {
            auto x = ec.path.at(k-1),
                 y = ec.path.at(k);
            auto e = gr->get_edge({x, y});
            for (uint64_t fr : e->frames) corr[fr] ^= 1;
            log_pr += log(e->probability);
#ifdef MEMORY_DEBUG
            std::cout << "\t\tApplied correction between " << print_v(x) << " and "
                << print_v(y) << ":";
            for (uint64_t fr : e->frames) std::cout << " " << fr;
            std::cout << std::endl;
#endif
        }
    }
    return log_pr;
}

}   // qontra
