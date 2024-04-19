/*
 *  author: Suhas Vittal
 *  date:   30 March 2024
 * */

#include "qontra/decoder/mobius.h"
#include "qontra/graph/decoding_graph/unified_lattice.h"

#include <PerfectMatching.h>

namespace qontra {

using namespace graph;
using namespace decoding;

MobiusDecoder::MobiusDecoder(const DetailedStimCircuit& circuit)
    :RestrictionDecoder(circuit),
    ufl_graph(),
    ufl_to_orig_map()
{
    ufl_graph = std::make_unique<DecodingGraph>(decoding_graph->make_unified_lattice(ufl_to_orig_map));
    ufl_graph->init_distances_for();
}

void
MobiusDecoder::load_syndrome(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome, int, int, bool) {
    RestrictionDecoder::load_syndrome(syndrome);
    // Update ufl_detectors.
    ufl_detectors.clear();
    for (uint64_t d : detectors) {
        if (d == get_color_boundary_index(COLOR_ANY)) {
            continue;
        }
        sptr<vertex_t> v = decoding_graph->get_vertex(d);
        for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
            for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
                if (v->color == c1 || v->color == c2) {
                    uint64_t _d = unified_lattice_id(v->id, c1, c2);
                    ufl_detectors.push_back(_d);
                }
            }
        }
    }
    // Activate detectors in ufl_graph.
    ufl_graph->activate_detectors(ufl_detectors, flags);
}

std::vector<assign_t>
MobiusDecoder::compute_matchings(stim::simd_bits_range_ref<SIMD_WIDTH>) {
    return compute_matching_on_unified_lattice();
}

std::vector<assign_t>
MobiusDecoder::compute_matching_on_unified_lattice() {
    const size_t n = ufl_detectors.size();
    const size_t m = n*(n+1)/2;

    PerfectMatching pm(n, m);
    pm.options.verbose = false;
    // There are no boundaries in the unified lattice, as there are an even number of detectors.
#ifdef MEMORY_DEBUG
    std::cout << "UFL Weights:" << std::endl;
#endif
    for (size_t i = 0; i < n; i++) {
        uint64_t di = ufl_detectors.at(i);
#ifdef MEMORY_DEBUG
        std::cout << "\t" << print_ufl(di) << ":";
#endif
        for (size_t j = i+1; j < n; j++) {
            uint64_t dj = ufl_detectors.at(j);
            error_chain_t ec = ufl_graph->get_error_chain(di, dj);
            uint32_t iw = static_cast<uint32_t>(1000 * ec.weight);
            pm.AddEdge(i, j, iw);
#ifdef MEMORY_DEBUG
            std::cout << " " << print_ufl(dj) << "[ P = " << ec.probability << " ]";
#endif
        }
#ifdef MEMORY_DEBUG
        std::cout << std::endl;
#endif
    }
    pm.Solve();
    // Return assignments. Note that these must be translated back to the original vertices.
    std::vector<assign_t> assign_arr;
    for (size_t i = 0; i < n; i++) {
        size_t j = pm.GetMatch(i);
        uint64_t di = ufl_detectors.at(i),
                 dj = ufl_detectors.at(j);
        if (di > dj) continue;
        sptr<vertex_t> x = ufl_graph->get_vertex(di),
                        y = ufl_graph->get_vertex(dj);
        read_ufl_error_chain(assign_arr, x, y);
    }
    for (assign_t& m : assign_arr) {
        identify_flag_edges_in_path(m);
    }
#ifdef MEMORY_DEBUG
    for (const assign_t& m : assign_arr) {
        error_chain_t ec = decoding_graph->get_error_chain(m.v, m.w, m.c1, m.c2);
        error_chain_t _ec = decoding_graph->get_error_chain(m.v, m.w, m.c1, m.c2, true);
        if (ec.path.front() != m.v) std::reverse(ec.path.begin(), ec.path.end());
        if (ec.path != m.path) {
            std::cout << "Path mismatch! mobius path: [";
            for (sptr<vertex_t> v : m.path) std::cout << " " << print_v(v);
            std::cout << " ], flagged path: [";
            for (sptr<vertex_t> v : ec.path) std::cout << " " << print_v(v);
            std::cout << " ], unflagged path: [";
            for (sptr<vertex_t> v : _ec.path) std::cout << " " << print_v(v);
            std::cout << " ]" << std::endl;
        }
    }
#endif
    return assign_arr;
}

void
MobiusDecoder::read_ufl_error_chain(std::vector<assign_t>& arr, sptr<vertex_t> src, sptr<vertex_t> dst) {
    assign_t curr;
    curr.v = ufl_to_orig_map.at(src);
    curr.w = nullptr;
    curr.path = { ufl_to_orig_map.at(src) };

    get_unified_lattice_bits_from_id(src->id, curr.c1, curr.c2);

    error_chain_t ec = ufl_graph->get_error_chain(src, dst);
    if (ec.path.empty()) return;
    if (ec.path.front() != src) std::reverse(ec.path.begin(), ec.path.end());

#ifdef MEMORY_DEBUG
    error_chain_t _ec = ufl_graph->get_error_chain(src, dst, COLOR_ANY, COLOR_ANY, true);
    std::cout << "UFL path between " << print_ufl(src) << " and " << print_ufl(dst) << ": flagged = [";
    for (sptr<vertex_t> v : ec.path) std::cout << " " << print_ufl(v);
    std::cout << " ], unflagged = [";
    for (sptr<vertex_t> v : _ec.path) std::cout << " " << print_ufl(v);
    std::cout << " ]" << std::endl;
#endif

    for (size_t i = 1; i < ec.path.size(); i++) {
        sptr<vertex_t> x = ec.path[i-1],
                        y = ec.path[i];
        sptr<vertex_t> v = ufl_to_orig_map.at(x),
                        w = ufl_to_orig_map.at(y);
        // If (x, y) is a crease edge, then we need to match to a boundary.
        int xc1, xc2, yc1, yc2;
        uint64_t lx = get_unified_lattice_bits_from_id(x->id, xc1, xc2),
                 ly = get_unified_lattice_bits_from_id(y->id, yc1, yc2);
        if (lx != ly) { // Crease edge:
            // Check which boundary it is by finding the common color between xc1, xc2, yc1, and yc2.
            int bcolor = (xc1 == yc1 || xc1 == yc2) ? xc1 : xc2;
            sptr<vertex_t> bb = decoding_graph->get_boundary_vertex(bcolor);
            curr.w = bb;
            curr.path.push_back(bb);
            arr.push_back(curr);
            // Make a new assignment now.
            curr.v = bb;
            curr.w = w;
            curr.c1 = yc1;
            curr.c2 = yc2;
            curr.path = { bb, w };
        } else {
            curr.w = w;
            curr.path.push_back(w);
        }
    }
    if (curr.v != nullptr && curr.w != nullptr) {
        arr.push_back(curr);
    }
}


}   // qontra
