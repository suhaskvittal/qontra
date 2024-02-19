/*
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#define MEMORY_DEBUG

#include "qontra/decoder/restriction.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#include <initializer_list>

const int COLOR_RED = 0;

namespace qontra {

using namespace graph;
using namespace decoding;

template <class SETLIKE> inline size_t
locally_matches(SETLIKE s1, const std::set<vpair_t>& s2, sptr<vertex_t> v) {
    size_t m = 0;
    for (const vpair_t& e : s2) {
        if (!s1.count(e)) return 0;
        else              m++;
    }
    return m;
}

inline bool
update_best_boundary(
        int intersection,
        fp_t& best_p,
        std::set<vpair_t>& best_boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> best_corr,
        fp_t& p,
        std::set<vpair_t>& boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr)
{
    if (intersection > 0 && p > best_p) {
        best_p = p;
        best_boundary = boundary;
        best_corr.clear();
        best_corr |= corr;
        return true;
    }
    return false;
}

void
remove_widowed_edges(std::map<vpair_t, size_t>& incidence_map) {
    std::map<sptr<vertex_t>, size_t> vertex_inc_map;
    for (auto& p : incidence_map) {
        const vpair_t& e = p.first;
        sptr<vertex_t> v1 = e.first,
                       v2 = e.second;
        vertex_inc_map[v1]++;
        vertex_inc_map[v2]++;
    }
    // Remove any pairs of vertices where both endpoints only have a single
    // incidence.
    for (auto it = incidence_map.begin(); it != incidence_map.end(); ) {
        vpair_t e = it->first;
        sptr<vertex_t> v1 = e.first,
                       v2 = e.second;
        if (vertex_inc_map.at(v1) == 1 && vertex_inc_map.at(v2) == 1) {
            it = incidence_map.erase(it);
        } else {
            it++;
        }
    }
}

Decoder::result_t
RestrictionDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();
    stim::simd_bits<SIMD_WIDTH> corr(n_obs);

    load_syndrome(syndrome);
    auto _detectors = std::move(detectors);
    auto _flags = std::move(flags);
    if (_detectors.empty()) return { 0.0, corr };
    // Compute the MWPM for each restricted lattice.
    std::vector<assign_t> matchings;
    for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
        for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
            load_syndrome(syndrome, c1, c2);
            std::vector<Decoder::assign_t> _matchings = compute_matching(c1, c2, true);
#ifdef MEMORY_DEBUG
            std::cout << "\tMatchings:" << std::endl;
#endif
            for (Decoder::assign_t x : _matchings) {
#ifdef MEMORY_DEBUG
                std::cout << "\t\t" << std::get<0>(x) << " <---> " << std::get<1>(x) << std::endl;
#endif
                matchings.push_back(cast_assign(x, c1, c2));
            }
        }
    }
    // Compute connected components.
    std::vector<component_t> components = compute_connected_components(matchings);
    // Identify all edges inside a connected component (and which component it is), and
    // do the same for all edges outside a connected component. Furthermore, maintain a
    // counter for the number of times an edge appears in a connected component.
    std::map<vpair_t, size_t> in_cc_map, not_cc_map;
    std::vector<std::pair<std::set<vpair_t>, int>> component_edge_sets;
    // Track all matches that are in a component.
    std::set<assign_t> in_cc_assignments;
#ifdef MEMORY_DEBUG
    std::cout << "In components:" << std::endl;
#endif
    for (const component_t& cc : components) {
        int color = cc.color;
        std::set<vpair_t> edge_set;
        for (const assign_t& m : cc.assignments) {
            in_cc_assignments.insert(m);
            sptr<vertex_t> v = decoding_graph->get_vertex(std::get<0>(m)),
                           w = decoding_graph->get_vertex(std::get<1>(m));
#ifdef MEMORY_DEBUG
            std::cout << "\t" << print_v(v) << " <---> " << print_v(w) << " on RL(" << std::get<2>(m) << "," << std::get<3>(m) << ")" << std::endl;
#endif
            std::set<vpair_t> tmp =
                insert_error_chain_into(in_cc_map, v, w, color, std::get<2>(m), std::get<3>(m));
            vtils::insert_range(edge_set, tmp);
        }
        component_edge_sets.emplace_back(edge_set, color);
    }
    // Do not_cc now.
#ifdef MEMORY_DEBUG
    std::cout << "Not in any component:" << std::endl;
#endif
    for (const assign_t& m : matchings) {
        if (in_cc_assignments.count(m)) continue;
        sptr<vertex_t> v = decoding_graph->get_vertex(std::get<0>(m)),
                       w = decoding_graph->get_vertex(std::get<1>(m));
#ifdef MEMORY_DEBUG
            std::cout << "\t" << print_v(v) << " <---> " << print_v(w) << " on RL(" << std::get<2>(m) << "," << std::get<3>(m) << ")" << std::endl;
#endif
        insert_error_chain_into(not_cc_map, v, w, COLOR_RED, std::get<2>(m), std::get<3>(m));
    }

    if (in_cc_map.size() && not_cc_map.size()) {
        return { 0, corr };
    }
    // Here, we will iterate multiple times and try to match faces as much as possible.
    size_t tries = 0;
r_compute_correction:
    std::set<sptr<vertex_t>> not_cc_incident,
                             in_cc_incident;
    insert_incident_vertices(not_cc_incident, not_cc_map, COLOR_RED);
    for (auto& es : component_edge_sets) {
        // Remove any deleted edges from the edge set.
        for (auto it = es.first.begin(); it != es.first.end(); ) {
            if (!in_cc_map.count(*it))  it = es.first.erase(it);
            else                        it++;
        }
        insert_incident_vertices(in_cc_incident, es.first, es.second);
    }
    std::set<sptr<vertex_t>> all_incident(not_cc_incident);
    vtils::insert_range(all_incident, in_cc_incident);
#ifdef MEMORY_DEBUG
    std::cout << "Incident vertices: "
        << not_cc_incident.size() << ", "
        << in_cc_incident.size() << ", "
        << all_incident.size() << std::endl;
#endif

    for (sptr<vertex_t> v : all_incident) {
        std::set<sptr<hyperedge_t>> faces = get_faces(v);
        const size_t nf = faces.size();
        const uint64_t enf = 1L << nf;
#ifdef MEMORY_DEBUG
        std::cout << "Faces (nf = " << nf << "):" << std::endl;
        for (auto f : faces) {
            std::cout << "\t<";
            for (size_t j = 0; j < f->get_order(); j++) {
                std::cout << " " << print_v(f->get<vertex_t>(j));
            }
            std::cout << " >, corr =";
            for (uint64_t fr : f->frames) {
                std::cout << " " << fr;
            }
            std::cout << std::endl;
        }
#endif
        // Track intersections with connected components and outside of
        // connected components.
        std::set<vpair_t> best_cc_boundary,
                          best_no_cc_boundary;
        stim::simd_bits<SIMD_WIDTH> best_cc_corr(n_obs),
                                    best_no_cc_corr(n_obs);
        fp_t best_prob_cc = 0.0,
             best_prob_no_cc = 0.0;
        for (uint64_t i = 0; i < enf; i++) {
            std::set<vpair_t> boundary;
            stim::simd_bits<SIMD_WIDTH> local_corr(n_obs);
            fp_t pr = 1.0;

            uint64_t ii = i;
            for (auto it = faces.begin(); it != faces.end() && ii; it++) {
                if (ii & 1) {
                    intersect_with_boundary(boundary, local_corr, pr, *it, v);
                }
                ii >>= 1;
            }
            // Now, we need to check how much boundary intersects with connected components
            // or anything not in the connected components. Note that this intersection is
            // in the neighborhood of v.
            size_t int_in_cc = locally_matches(in_cc_map, boundary, v),
                   int_not_cc = locally_matches(not_cc_map, boundary, v);
            update_best_boundary(
                int_in_cc, best_prob_cc, best_cc_boundary, best_cc_corr, pr, boundary, local_corr);
            update_best_boundary(
                int_not_cc, best_prob_no_cc, best_no_cc_boundary, best_no_cc_corr, pr, boundary, local_corr);
        }
        if (best_prob_cc > best_prob_no_cc) {
            corr ^= best_cc_corr;
        } else {
            corr ^= best_no_cc_corr;
        }
    }
    // Remove any widowed edges, as these can cause the decoder to loop infinitely.
    remove_widowed_edges(in_cc_map);
    remove_widowed_edges(not_cc_map);
    if (in_cc_map.size() > 1 && not_cc_map.size() > 1) {
        if (tries < 100) {
            tries++;
            goto r_compute_correction;
        } else {
            // Throw an error here.
            std::cerr << "[ RestrictionDecoder ] Failed to compute correction for syndrome: D[";
            for (uint64_t d : _detectors) std::cerr << " " << d;
            std::cerr << " ] and F[";
            for (uint64_t f : _flags) std::cerr << " " << f;
            std::cerr << " ]" << std::endl
                    << "\tEdges remaining in CC: " << in_cc_map.size() << std::endl
                    << "\tEdges remaining out of CC: " << not_cc_map.size() << std::endl;
            exit(1);
        }
    }
    // Otherwise, we are done.
    return { 0.0, corr };
}

std::vector<RestrictionDecoder::component_t>
RestrictionDecoder::compute_connected_components(
        const std::vector<RestrictionDecoder::assign_t>& assignments)
{
    // For each connected component, we know the following:
    //  (1) A boundary may be connected to multiple vertices (direct matches
    //      to the boundary, or matches that go through a boundary).
    //  (2) Non-boundary vertices are only connected to R-1 vertices (at most),
    //      where R is the number of restricted lattices.
    //
    // We can identify connected components via DFS from the red boundary.
    //
    // Create a graph so we can launch this DFS.
    struct e_t : base::edge_t {
        int c1;
        int c2;
    };
    auto cgr = std::make_unique<Graph<vertex_t, e_t>>();
    // Add all assignments to the graph.
    for (const assign_t& m : assignments) {
        uint64_t dv = std::get<0>(m),
                 dw = std::get<1>(m);
        sptr<vertex_t> v = decoding_graph->get_vertex(dv),
                       w = decoding_graph->get_vertex(dw);
        if (!cgr->contains(v)) cgr->add_vertex(v);
        if (!cgr->contains(w)) cgr->add_vertex(w);
        if (cgr->contains(v, w)) continue;
        sptr<e_t> e = cgr->make_and_add_edge(v, w);
        e->c1 = std::get<2>(m);
        e->c2 = std::get<3>(m);
    }
    // If the red boundary is not in the graph, then just exit as there will be
    // no connected components.
    sptr<vertex_t> vrb = decoding_graph->get_boundary_vertex(COLOR_RED);
    if (!cgr->contains(vrb)) return {};
    // Now, we will perform a DFS from vrb. This is a weird DFS in that:
    // vertices are NOT marked as visited. Instead, traversal will mark
    // edges on the graph. Once a boundary is hit, the marked edges will
    // be deleted. Any branches stemming from the deleted edges will also
    // be ignored.
    //
    // prev_edge_map keeps track of the previous edge in the traversal.
    // incoming_edge_map tracks what edge visited a vertex most recently in
    // the DFS.
    std::vector<component_t> components;

    std::map<sptr<e_t>, sptr<e_t>> prev_edge_map;
    std::map<sptr<vertex_t>, sptr<e_t>> incoming_edge_map;

    std::vector<sptr<vertex_t>> dfs{vrb};
    incoming_edge_map[vrb] = nullptr;
    while (dfs.size()) {
        sptr<vertex_t> v = dfs.back();
        dfs.pop_back();
#ifdef MEMORY_DEBUG
        std::cout << "traversing to " << print_v(v) << std::endl;
#endif
        // Make sure we are fine to continue.
        sptr<e_t> e = incoming_edge_map[v];
        if (e != nullptr && prev_edge_map[e] != nullptr && !cgr->contains(prev_edge_map[e])) {
            continue;
        }
        // If v is a boundary vertex, stop and identify the connected component.
        if (v->is_boundary_vertex && e != nullptr) {
            std::vector<assign_t> m_in_cc;
            sptr<e_t> curr = e;
            while (curr != nullptr) {
                sptr<vertex_t> x = e->get_source<vertex_t>(),
                                y = e->get_target<vertex_t>();
                m_in_cc.emplace_back(x->id, y->id, e->c1, e->c2);
                // Delete the edge from cgr.
                cgr->delete_edge(curr);
                curr = prev_edge_map[curr];
            }
            std::initializer_list<sptr<vertex_t>> blist{v, vrb};
            int cc_color =
                decoding_graph->get_complementary_boundaries_to(blist)[0]->color;
            components.push_back({m_in_cc, cc_color});
            continue;
        }
        // Visit neighbors (not in incoming edge).
#ifdef MEMORY_DEBUG
        std::cout << "\tneighbors:";
#endif
        for (sptr<vertex_t> w : cgr->get_neighbors(v)) {
            if (e != nullptr &&
                (e->get_source<vertex_t>() == w || e->get_target<vertex_t>() == w)) 
            {
                continue;
            }
            if (incoming_edge_map[w] != nullptr && cgr->contains(incoming_edge_map[w])) {
                continue;
            }
            dfs.push_back(w);
            sptr<e_t> _e = cgr->get_edge(v, w);
            incoming_edge_map[w] = _e;
            prev_edge_map[_e] = e;
#ifdef MEMORY_DEBUG
            std::cout << " " << print_v(w);
#endif
        }
#ifdef MEMORY_DEBUG
        std::cout << std::endl;
#endif
    }
    return components;
}

std::set<vpair_t>
RestrictionDecoder::insert_error_chain_into(
        std::map<vpair_t, size_t>& incidence_map,
        sptr<vertex_t> src,
        sptr<vertex_t> dst,
        int component_color,
        int c1,
        int c2)
{
    error_chain_t ec = decoding_graph->get(c1, c2, src, dst);
    for (size_t i = 1; i < ec.path.size(); i++) {
        sptr<vertex_t> v = ec.path.at(i-1),
                       w = ec.path.at(i);
        if (v->is_boundary_vertex && w->is_boundary_vertex) continue;
        if (v->color != component_color && w->color != component_color) continue;
#ifdef MEMORY_DEBUG
        std::cout << "\t\t" << print_v(v) << " , " << print_v(w);
#endif
        // Flatten v and w.
        sptr<vertex_t> fv = v->get_base(),
                       fw = w->get_base();
        std::cout << " --> " << print_v(fv) << " , " << print_v(fw) << std::endl;
        if (fv == fw) continue;
        // Check if fv and fw share an edge.
        if (!decoding_graph->share_hyperedge({fv, fw})) {
            // Then perhaps there is a CX edge. To account for this, add the path between
            // fv and fw.
            insert_error_chain_into(incidence_map, fv, fw, component_color, c1, c2);
        } else {
            // Make sure that their colors are not equal. Otherwise, get a common neighbor
            // of fv and fw which has a different color (call this fu). 
            // Add (fv, fu) and (fu, fw) instead.
            if (fv->color == fw->color) {
                sptr<vertex_t> fu = nullptr;
                for (sptr<vertex_t> x : decoding_graph->get_common_neighbors({v, w})) {
                    sptr<vertex_t> fx = x->get_base();
                    if (fv->color == fx->color) continue;
                    fu = fx;
                    break;
                }
                if (fu == nullptr) {
                    // This should not happen: exit.
                    std::cerr << "[ insert_error_chain_into ] could not find common vertex between "
                        << print_v(v) << "(" << print_v(fv) << ") and "
                        << print_v(w) << "(" << print_v(fw) << ")" << std::endl;
                    exit(1);
                }
                vpair_t e1 = make_vpair(fv, fu),
                        e2 = make_vpair(fu, fw);
                incidence_map[e1]++;
                incidence_map[e2]++;
            } else {
                vpair_t e = make_vpair(fv, fw);
                incidence_map[e]++;
            }
        }
    }
    std::set<vpair_t> edge_set;
    for (auto& p : incidence_map) edge_set.insert(p.first);
    return edge_set;
}

std::set<sptr<hyperedge_t>>
RestrictionDecoder::get_faces(sptr<vertex_t> v) {
    std::set<sptr<hyperedge_t>> faces;
    // Here, we will flatten the faces themselves. We require that
    // each face has different colors.
    for (sptr<hyperedge_t> e : decoding_graph->get_common_hyperedges({v})) {
        std::vector<sptr<vertex_t>> flat_vlist;
        std::set<int> colors_in_face;
        bool do_not_add = false;
        for (size_t i = 0; i < e->get_order(); i++) {
            sptr<vertex_t> x = e->get<vertex_t>(i);
            if (colors_in_face.count(x->color)
                || std::find(flat_vlist.begin(), flat_vlist.end(), x) != flat_vlist.end()) 
            {
                do_not_add = true;
                break;
            }
            flat_vlist.push_back(x->get_base());
            colors_in_face.insert(x->color);
        }
        if (do_not_add) continue;
        sptr<hyperedge_t> fe = decoding_graph->get_edge(flat_vlist);
        if (fe != nullptr) {
            faces.insert(fe);
        }
    }
    return faces;
}

}   // qontra
